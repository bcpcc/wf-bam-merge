#!/bin/bash -euo pipefail
mkdir -p 'out/unaligned_bams'
    echo '{
    "barcode": null,
    "type": "test_sample",
    "run_ids": [
        "267479e8-b237-4c66-8246-1cbfe2876eeb"
    ],
    "basecall_models": [
        "dna_r10.4.1_e8.2_400bps_sup@v4.3.0"
    ],
    "alias": "unaligned_bams",
    "src_xam": null,
    "src_xai": null,
    "is_unaligned": true,
    "n_primary": 0,
    "n_unmapped": 5000
}' > metamap.json
    mv metamap.json reads/reads.bam stats/bamstats_results index/reads.bam.bai 'out/unaligned_bams'
