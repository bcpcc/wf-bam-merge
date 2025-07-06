#!/bin/bash -euo pipefail
samtools cat -b <(find input_bams -name 'reads*.bam' | sort)     | samtools sort - -@ 2 --write-index -o reads.bam##idx##reads.bam.bai
