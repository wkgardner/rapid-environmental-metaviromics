#!/usr/bin/env python3
"""Upload a REMV results bundle to a LabKey experiment folder over WebDAV.

The analytical bundle is small: the per-sample TSV tables, the MultiQC report and the run
metadata, about 6 MB in total for a 48-sample run. EsViritu's BAM files are excluded unless
asked for, because they are tens of gigabytes and WebDAV is a poor way to move them.

This follows the same WebDAV mechanics as NVD's bin/webdav_CLIent.py: HTTP Basic auth with
the literal username "apikey", and a plain PUT per file. As NVD does, it never calls MKCOL:
LabKey's WebDAV creates the parent collections itself as the first file is written into
them, and many accounts that can write files are refused an explicit MKCOL.

It differs from NVD in one deliberate way: the API key is read from a file, never from an
environment variable, so it cannot leak into a process listing, a shell history or a
Nextflow task environment.

Nothing is uploaded unless --upload is passed. The default is a dry run that prints exactly
what would be sent.

    # see what would happen, touches nothing
    pixi run python scripts/upload_results_to_labkey.py \
        --results /staging/groups/oconnor_group/rapid-ww-metagenomics/32824/results

    # actually upload
    pixi run python scripts/upload_results_to_labkey.py \
        --results /staging/groups/oconnor_group/rapid-ww-metagenomics/32824/results \
        --api-key-file ~/.labkey/dholk_api_key \
        --upload
"""

from __future__ import annotations

import argparse
import base64
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

# The experiment's WebDAV root on DholK. The browsable page for the same container is
# https://dholk.primate.wisc.edu/dho/experiments/32824/project-begin.view
DEFAULT_WEBDAV_BASE = "https://dholk.primate.wisc.edu/_webdav/dho/experiments/32824/@files/"
DEFAULT_REMOTE_FOLDER = "REMv-results"

# LabKey expects the literal string "apikey" as the username when the password is an API key.
WEBDAV_USERNAME = "apikey"

# The bundle, as (glob relative to the results folder, what it is). Order is the order of
# upload, so the cheap metadata lands before the bulkier per-sample files.
#
# This covers the pipeline's own tables AND EsViritu's native per-sample output. EsViritu
# writes nine files per sample beside its temp folder; uploading only the call table loses
# the coverage windows, the consensus sequences and the run's own parameters.
BUNDLE_PATTERNS: tuple[tuple[str, str], ...] = (
    ("05_run_summary/*.run_metadata.yaml", "run metadata"),
    ("04_qc_report/*.multiqc_report.html", "MultiQC report"),
    ("01_preprocessing/preprocessing_read_counts/*.tsv", "preprocessing read counts"),
    ("02_genome_kmer_based_calling/sylph_results/*.tsv", "sylph profiles"),
    ("05_run_summary/sample_status/*.tsv", "sample status"),
    # EsViritu's own per-sample output: tables, consensus, log, params.
    ("03_genome_mapping_based_calling/esviritu_results/*/*.tsv", "EsViritu tables"),
    ("03_genome_mapping_based_calling/esviritu_results/*/*.fasta", "EsViritu consensus FASTA"),
    ("03_genome_mapping_based_calling/esviritu_results/*/*.yaml", "EsViritu params/readstats"),
    ("03_genome_mapping_based_calling/esviritu_results/*/*.json", "EsViritu fastp JSON"),
    ("03_genome_mapping_based_calling/esviritu_results/*/*.log", "EsViritu log"),
    # The round-3 alignments. These are the reviewer's route back to raw reads, so they ship by
    # default; about 8.7 GB for 48 samples. The dry run prints the total before anything moves.
    ("03_genome_mapping_based_calling/esviritu_final_bams/*.bam", "EsViritu final BAM"),
    ("03_genome_mapping_based_calling/esviritu_final_bams/*.bam.bai", "EsViritu final BAM index"),
    # The reference set round 3 mapped against. Small, and without it the BAM supports
    # idxstats/view/fastq but not pileup, variant calling or recomputing identity.
    ("03_genome_mapping_based_calling/esviritu_final_bams/*.third.reference.fa",
     "EsViritu round-3 reference"),
)

# --include-clustering. The read-sharing and clustering tables EsViritu keeps alongside its
# earlier alignments, which are what the 33% read-sharing cutoff actually operates on. About
# 320 MB for 48 samples, nearly all of it the *_read_comp_compare.tsv pair.
CLUSTERING_PATTERNS: tuple[tuple[str, str], ...] = (
    ("03_genome_mapping_based_calling/esviritu_intermediates/*/*.tsv",
     "EsViritu clustering tables"),
    ("03_genome_mapping_based_calling/esviritu_intermediates/*/*.fasta",
     "EsViritu clustered reference FASTA"),
)

# --include-bams. The round-1 and round-2 alignments, ~49 GB; WebDAV is a poor way to move
# these. Only needed to reproduce the clustering dropout analysis.
BAM_PATTERNS: tuple[tuple[str, str], ...] = (
    ("03_genome_mapping_based_calling/esviritu_intermediates/*/*.bam", "EsViritu round 1-2 BAM"),
    ("03_genome_mapping_based_calling/esviritu_intermediates/*/*.bam.bai",
     "EsViritu round 1-2 BAM index"),
)


def read_api_key(api_key_file: Path) -> str:
    """Read the API key from its own file.

    The key is never echoed, never put in the environment and never written to a log. The
    file should be readable only by its owner; this refuses to continue if it is not.
    """
    resolved = api_key_file.expanduser()
    if not resolved.is_file():
        raise SystemExit(f"API key file not found: {resolved}")

    mode = resolved.stat().st_mode
    if mode & 0o077:
        raise SystemExit(
            f"API key file is readable by others: {resolved}\n"
            f"Fix it with: chmod 600 {resolved}"
        )

    # Take the first line that is neither blank nor a comment, so a file that documents
    # itself still works.
    candidate = ""
    for line in resolved.read_text().splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            candidate = stripped
            break

    # Accept an env-style line such as LABKEY_API_KEY=apikey|abc... as well as a bare key.
    # Sending the whole "NAME=value" string as the password is a silent HTTP 401, which is
    # indistinguishable from an expired key unless it is handled here.
    if "=" in candidate and not candidate.startswith("apikey|"):
        candidate = candidate.split("=", 1)[1].strip()

    # Shell-style quoting survives copy and paste surprisingly often.
    if len(candidate) >= 2 and candidate[0] == candidate[-1] and candidate[0] in "\"'":
        candidate = candidate[1:-1].strip()

    if not candidate:
        raise SystemExit(f"API key file holds no usable key: {resolved}")
    return candidate


def build_auth_header(api_key: str) -> str:
    """Build the HTTP Basic auth header value for LabKey's apikey scheme."""
    credentials = f"{WEBDAV_USERNAME}:{api_key}"
    return "Basic " + base64.b64encode(credentials.encode()).decode()


def collect_bundle_files(
    results_directory: Path,
    include_clustering: bool,
    include_bams: bool,
) -> list[tuple[Path, str, str]]:
    """Find every file in the bundle.

    Returns (local path, remote path relative to the bundle folder, description). The remote
    path mirrors the published results layout, so the LabKey folder reads the same way the
    results directory does.
    """
    patterns = BUNDLE_PATTERNS
    if include_clustering:
        patterns += CLUSTERING_PATTERNS
    if include_bams:
        patterns += BAM_PATTERNS

    # A path can match more than one pattern, so remember what has been claimed and keep the
    # first description. Without this a file would be uploaded, and counted, twice.
    found: list[tuple[Path, str, str]] = []
    already_claimed: set[Path] = set()

    for pattern, description in patterns:
        for local_path in sorted(results_directory.glob(pattern)):
            if local_path.is_file() and local_path not in already_claimed:
                already_claimed.add(local_path)
                relative_path = local_path.relative_to(results_directory).as_posix()
                found.append((local_path, relative_path, description))

    return found


def webdav_request(
    url: str,
    method: str,
    auth_header: str,
    payload: bytes | None = None,
) -> int:
    """Send one WebDAV request and return its HTTP status code.

    HTTP errors are returned as their status code rather than raised, because the useful
    WebDAV responses include several that urllib treats as errors: 405 means a collection
    already exists, 404 means it does not yet.
    """
    request = urllib.request.Request(url, data=payload, method=method)
    request.add_header("Authorization", auth_header)
    try:
        with urllib.request.urlopen(request) as response:
            return response.getcode()
    except urllib.error.HTTPError as error:
        return error.code
    except urllib.error.URLError as error:
        raise SystemExit(f"Could not reach {url}: {error.reason}") from error


def create_collection(url: str, auth_header: str) -> bool:
    """Try to create a WebDAV collection. Only used as a fallback; see upload_one_file.

    405 (method not allowed) is LabKey's answer for a collection that already exists, so it
    counts as success. Many LabKey accounts can write files but are not allowed to MKCOL at
    all, which is why this is never on the main path.
    """
    directory_url = url if url.endswith("/") else url + "/"
    return webdav_request(directory_url, "MKCOL", auth_header) in (201, 204, 405)


def upload_one_file(local_path: Path, remote_url: str, auth_header: str) -> tuple[bool, int]:
    """PUT a single file, returning (accepted, final status code).

    LabKey's WebDAV creates missing parent collections on PUT, so no MKCOL is needed and
    none is sent: NVD's client uploads several levels deep without ever calling mkdir, and
    an explicit MKCOL is refused for accounts that can nonetheless write files.

    A 409 is the one case where the server really does want the parent first, so that alone
    falls back to creating parents and retrying once.
    """
    payload = local_path.read_bytes()
    status = webdav_request(remote_url, "PUT", auth_header, payload=payload)
    if status in (200, 201, 204):
        return True, status

    if status == 409:
        parent_url = remote_url.rsplit("/", 1)[0] + "/"
        if create_collection(parent_url, auth_header):
            status = webdav_request(remote_url, "PUT", auth_header, payload=payload)
            return status in (200, 201, 204), status

    return False, status


def format_size(total_bytes: int) -> str:
    """Human-readable byte count."""
    size = float(total_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Upload a REMV results bundle to LabKey over WebDAV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--results", type=Path, required=True,
                        help="The run's published results directory.")
    parser.add_argument("--webdav-base", default=DEFAULT_WEBDAV_BASE,
                        help=f"Experiment WebDAV root. Default: {DEFAULT_WEBDAV_BASE}")
    parser.add_argument("--remote-folder", default=DEFAULT_REMOTE_FOLDER,
                        help=f"Folder to create under the root. Default: {DEFAULT_REMOTE_FOLDER}")
    parser.add_argument("--api-key-file", type=Path,
                        help="File holding the LabKey API key. Required with --upload.")
    parser.add_argument("--include-clustering", action="store_true",
                        help="Also upload EsViritu's read-sharing and clustering tables from "
                             "esviritu_intermediates/ (about 320 MB for 48 samples).")
    parser.add_argument("--include-bams", action="store_true",
                        help="Also upload EsViritu's round-1 and round-2 alignments from "
                             "esviritu_intermediates/ (~49 GB). The final round-3 BAMs are in "
                             "the default bundle already; this is only for reproducing the "
                             "clustering dropout analysis.")
    parser.add_argument("--upload", action="store_true",
                        help="Actually upload. Without this, the script only reports what it would do.")
    parser.add_argument("--check-auth", action="store_true",
                        help="Test the API key with one harmless GET and exit. Writes nothing.")
    args = parser.parse_args()

    if args.check_auth:
        if not args.api_key_file:
            raise SystemExit("--check-auth requires --api-key-file")
        base = args.webdav_base if args.webdav_base.endswith("/") else args.webdav_base + "/"
        auth_header = build_auth_header(read_api_key(args.api_key_file))
        status = webdav_request(base, "GET", auth_header)
        print(f"GET {base} -> HTTP {status}")
        if status == 200:
            print("Key accepted and the container is readable.")
            return 0
        if status == 401:
            print("HTTP 401: the key itself was not accepted. It is wrong, expired, or the\n"
                  "file holds something other than the key. Note that NVD does not read a\n"
                  "key file at all: it uses the Nextflow secret LABKEY_API_KEY.")
            return 1
        if status == 403:
            print("HTTP 403: the key is valid but lacks permission on this container.")
            return 1
        print("Unexpected status; see above.")
        return 1

    results_directory = args.results.expanduser().resolve()
    if not results_directory.is_dir():
        raise SystemExit(f"Results directory not found: {results_directory}")

    bundle = collect_bundle_files(results_directory, args.include_clustering, args.include_bams)
    if not bundle:
        raise SystemExit(f"No bundle files found under {results_directory}")

    total_bytes = sum(local.stat().st_size for local, _, _ in bundle)
    base_url = args.webdav_base if args.webdav_base.endswith("/") else args.webdav_base + "/"
    bundle_root = f"{args.remote_folder.strip('/')}/"

    counts_by_description: dict[str, int] = {}
    for _, _, description in bundle:
        counts_by_description[description] = counts_by_description.get(description, 0) + 1

    print(f"Results:  {results_directory}")
    print(f"Target:   {urllib.parse.urljoin(base_url, bundle_root)}")
    print(f"Bundle:   {len(bundle)} files, {format_size(total_bytes)}")
    for description, count in counts_by_description.items():
        print(f"            {count:>5}  {description}")
    print()

    if not args.upload:
        print("Dry run. Nothing was sent. Re-run with --api-key-file and --upload to upload.")
        return 0

    if not args.api_key_file:
        raise SystemExit("--upload requires --api-key-file")

    auth_header = build_auth_header(read_api_key(args.api_key_file))

    # No folder is created up front. LabKey's WebDAV makes the parent collections as the
    # first file is PUT into them, and an explicit MKCOL is refused for accounts that can
    # still write files perfectly well.
    failures: list[tuple[Path, int]] = []
    for index, (local_path, relative_path, _) in enumerate(bundle, start=1):
        remote_url = urllib.parse.urljoin(base_url, f"{bundle_root}{relative_path}")
        accepted, status = upload_one_file(local_path, remote_url, auth_header)
        if accepted:
            print(f"  [{index}/{len(bundle)}] {relative_path}")
        else:
            print(f"  [{index}/{len(bundle)}] FAILED (HTTP {status}) {relative_path}",
                  file=sys.stderr)
            failures.append((local_path, status))
            if status in (401, 403):
                raise SystemExit(
                    f"Not authorised to write {remote_url} (HTTP {status}).\n"
                    "The API key may be expired, or it may lack insert permission on this "
                    "container. Nothing further was attempted."
                )

    print()
    if failures:
        print(f"{len(failures)} of {len(bundle)} files failed to upload.", file=sys.stderr)
        return 1

    print(f"Uploaded {len(bundle)} files to {urllib.parse.urljoin(base_url, bundle_root)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
