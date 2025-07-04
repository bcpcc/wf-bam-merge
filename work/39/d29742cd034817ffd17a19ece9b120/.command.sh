#!/bin/bash -euo pipefail
mkdir bamstats_results
bamstats reads.bam -s unaligned_bams -u         -f bamstats_results/bamstats.flagstat.tsv -t 2         -i bamstats_results/bamstats.runids.tsv         -l bamstats_results/bamstats.basecallers.tsv         --histograms histograms      > /dev/null
mv histograms/* bamstats_results/

# get n_seqs from flagstats - need to sum them up
awk 'NR==1{for (i=1; i<=NF; i++) {ix[$i] = i}} NR>1 {c+=$ix["total"]} END{print c}'         bamstats_results/bamstats.flagstat.tsv > bamstats_results/n_seqs
# get unique run IDs (we add `-F '\t'` as `awk` uses any stretch of whitespace
# as field delimiter otherwise and thus ignore empty columns)
awk -F '\t' '
    NR==1 {for (i=1; i<=NF; i++) {ix[$i] = i}}
    # only print run_id if present
    NR>1 && $ix["run_id"] != "" {print $ix["run_id"]}
' bamstats_results/bamstats.runids.tsv | sort | uniq > bamstats_results/run_ids
# get unique basecall models
awk -F '\t' '
    NR==1 {for (i=1; i<=NF; i++) {ix[$i] = i}}
    # only print run_id if present
    NR>1 && $ix["basecaller"] != "" {print $ix["basecaller"]}
' bamstats_results/bamstats.basecallers.tsv | sort | uniq > bamstats_results/basecallers
