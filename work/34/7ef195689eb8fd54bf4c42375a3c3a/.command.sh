#!/bin/bash -euo pipefail
# Output nextflow params object to JSON
    echo '{
    "help": false,
    "version": false,
    "bam": "test_data/unaligned_bams/",
    "out_dir": "test_results_unaligned",
    "sample": null,
    "sample_sheet": null,
    "client_fields": null,
    "alignment_status": "unaligned",
    "sort_before_merge": false,
    "aws_image_prefix": null,
    "aws_queue": null,
    "disable_ping": false,
    "monochrome_logs": false,
    "validate_params": true,
    "show_hidden_params": false,
    "schema_ignore_params": "show_hidden_params,validate_params,monochrome_logs,aws_queue,aws_image_prefix,wf,store_dir",
    "wf": {
        "bamstats": true,
        "keep_unaligned": true,
        "per_read_stats": false,
        "allow_multiple_basecall_models": false,
        "example_cmd": [
            "--bam '\''test_data/aligned_bams/'\'' --alignment_status aligned",
            "--bam '\''test_data/unaligned_bams/'\'' --alignment_status unaligned"
        ],
        "common_sha": "sha72f3517dd994984e0e2da0b97cb3f23f8540be4b",
        "agent": null,
        "epi2me_instance": null,
        "epi2me_user": null
    }
}' > params.json
