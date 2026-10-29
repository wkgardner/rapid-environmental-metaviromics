#!/usr/bin/env python3
"""Write <run_id>.run_metadata.yaml: what this run did, with which tools and settings.

Output tables carry tool-neutral names, so this file is where a reader finds which tool and
version ran each step, the reference filter's mode (deplete or extract) and index, which databases
the callers used, and every effective param.

--versions takes the versions.yml files the processes emit ('"WORKFLOW:...:PROCESS": {tool:
version}'); each process maps to one step. A step whose --skip_<step> param is true is recorded
as {skipped: true}; adapter trimming and dedup are switched off with NVD's trim_adapters / dedup
params instead.
"""

import argparse
import json
from pathlib import Path

import yaml

STEP_BY_PROCESS = {
    "REFERENCE_FILTER": "reference_filter",
    "REFERENCE_FILTER_SRA": "reference_filter",
    "BBDUK_TRIM_ADAPTERS": "adapter_trimming",
    "CLUMPIFY_DEDUP": "dedup",
    "BBDUK_FILTER_READS": "read_filtering",
    "SYLPH_SKETCH_READS": "genome_kmer_based_calling",
    "SYLPH_PROFILE": "genome_kmer_based_calling",
    "ESVIRITU": "genome_mapping_based_calling",
    "MULTIQC": "qc_report",
}
STEPS_IN_ORDER = ["reference_filter", "adapter_trimming", "dedup", "read_filtering",
                  "genome_kmer_based_calling", "genome_mapping_based_calling", "qc_report"]


def step_is_switched_off(step: str, params: dict) -> bool:
    nvd_on_off_param = {"adapter_trimming": "trim_adapters", "dedup": "dedup"}.get(step)
    if nvd_on_off_param is not None:
        return params.get(nvd_on_off_param) is False
    return bool(params.get(f"skip_{step}"))


def step_details(step: str, params: dict) -> dict:
    """What each step ran against, from the effective params."""
    if step == "reference_filter":
        return {"mode": params.get("reference_filter_mode"), "index": params.get("reference_filter_index"),
                "index_contents": params.get("reference_filter_index_description")}
    if step == "genome_kmer_based_calling":
        return {"database": params.get("sylph_database")}
    if step == "genome_mapping_based_calling":
        return {"database": params.get("esviritu_database")}
    return {}


def build_steps(version_files: list[Path], params: dict) -> dict:
    tool_versions_by_step: dict[str, dict] = {}
    for version_file in version_files:
        for process_name, tool_versions in (yaml.safe_load(version_file.read_text()) or {}).items():
            step = STEP_BY_PROCESS.get(process_name.split(":")[-1])
            if step is not None:
                tool, version = next(iter(tool_versions.items()))
                tool_versions_by_step[step] = {"tool": tool, "version": str(version)}

    steps = {}
    for step in STEPS_IN_ORDER:
        if step_is_switched_off(step, params):
            steps[step] = {"skipped": True}
        elif step in tool_versions_by_step:
            steps[step] = tool_versions_by_step[step] | step_details(step, params)
    return steps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--start-time", required=True)
    parser.add_argument("--params-json", required=True, type=Path, help="the run's effective params")
    parser.add_argument("--versions", nargs="*", default=[], type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    arguments = parser.parse_args()

    params = json.loads(arguments.params_json.read_text())
    metadata = {
        "run_id": arguments.run_id,
        "pipeline_version": arguments.pipeline_version,
        "start_time": arguments.start_time,
        "samplesheet": params.get("input"),
        "references": {key: params.get(key) for key in ("vpd_version", "vpd_mask_version", "taxdump_version")},
        "steps": build_steps(arguments.versions, params),
        "params": params,
    }
    arguments.out_dir.mkdir(parents=True, exist_ok=True)
    (arguments.out_dir / f"{arguments.run_id}.run_metadata.yaml").write_text(
        yaml.safe_dump(metadata, sort_keys=False))


if __name__ == "__main__":
    main()
