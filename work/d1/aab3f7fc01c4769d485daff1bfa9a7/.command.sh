#!/bin/bash -euo pipefail
echo "samtools,$(samtools --version | head -n1 | cut -d' ' -f2)" >> versions.txt
