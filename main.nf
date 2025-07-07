#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "wf_bam_merge"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt"
    cpus 1
    
    output:
    path "versions.txt"
    
    script:
    """
    echo "samtools,\$(samtools --version | head -n1 | cut -d' ' -f2)" >> versions.txt
    """
}
process SAMTOOLS_MERGE_ALIGNED {
    label "wf_bam_merge"
    publishDir "${params.out_dir}/merged_bams", mode: 'copy'
    
    input:
    tuple val(meta), path(bams), path(indices)
    
    output:
    tuple val(meta), path("*.merged.bam"), path("*.merged.bam.bai"), emit: merged_bam
    
    when:
    bams.size() > 1 && params.alignment_status == 'aligned'
    
    script:
    // Build output name with sample as default prefix
    def sample_name = params.sample ?: meta.alias ?: meta.id ?: "sample"
    def custom_prefix = params.output_prefix ?: sample_name
    def alignment_suffix = params.include_alignment_status ? ".aligned" : ""
    def date_suffix = params.include_date ? ".${new Date().format('yyyyMMdd')}" : ""
    def output_name = "${custom_prefix}${alignment_suffix}${date_suffix}.merged.bam"
    
    """
    echo "Original meta.alias: ${meta.alias}"
    echo "Using sample name: ${sample_name}"
    echo "Input BAMs (staged names): ${bams.join(' ')}"
    echo "Output will be: ${output_name}"
    
    # Merge aligned BAM files with samtools merge
    samtools merge \\
        -@ ${task.cpus} \\
        ${output_name} \\
        ${bams.join(' ')}
    
    # Index the merged BAM
    samtools index -@ ${task.cpus} ${output_name}
     echo "Successfully created: ${output_name}"
    ls -la *.merged.bam*
    """
}

process SAMTOOLS_CAT_UNALIGNED {
    label "wf_bam_merge"
    publishDir "${params.out_dir}/merged_bams", mode: 'copy'
    
    input:
    tuple val(meta), path(bams), path(indices)
    
    output:
    tuple val(meta), path("*.merged.bam"), path("*.merged.bam.bai"), emit: merged_bam
    
    when:
    bams.size() > 1 && params.alignment_status == 'unaligned'
    
    script:
    // Build output name with sample as default prefix
    def sample_name = params.sample ?: meta.alias ?: meta.id ?: "sample"
    def custom_prefix = params.output_prefix ?: sample_name
    def alignment_suffix = params.include_alignment_status ? ".unaligned" : ""
    def date_suffix = params.include_date ? ".${new Date().format('yyyyMMdd')}" : ""
    def output_name = "${custom_prefix}${alignment_suffix}${date_suffix}.merged.bam"
    
    """
    echo "Original meta.alias: ${meta.alias}"
    echo "Using sample name: ${sample_name}"
    echo "Input BAMs (staged names): ${bams.join(' ')}"
    echo "Output will be: ${output_name}"
    
    # Concatenate unaligned BAM files with samtools cat
    samtools cat \\
        -o ${output_name} \\
        ${bams.join(' ')}
    
    # Index the merged BAM
    samtools index -@ ${task.cpus} ${output_name}
     echo "Successfully created: ${output_name}"
    ls -la *.merged.bam*
    """
}

process SAMTOOLS_SORT {
    label "wf_bam_merge"
    
    input:
    tuple val(meta), path(bam), path(index)
    
    output:
    tuple val(meta), path("${meta.alias}.sorted.bam"), path("${meta.alias}.sorted.bam.bai"), emit: sorted_bam
    
    when:
    params.sort_before_merge && params.alignment_status == 'unaligned'
    
    script:
    def args = task.ext.args ?: ''
    def memory = task.memory ? "-m ${task.memory.toGiga()}G" : "-m 2G"
    
    """
    # Sort BAM with samtools sort
    samtools sort \\
        -@ ${task.cpus} \\
        ${memory} \\
        ${args} \\
        -o ${meta.alias}.sorted.bam \\
        ${bam}
    
    # Index the sorted BAM
    samtools index -@ ${task.cpus} ${meta.alias}.sorted.bam
    """
}

process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "wf-bam-merge*-report.html"
    
    input:
    tuple val(analysis_group), val(metadata), path(merged_bams, stageAs: "merged_*")
    path client_fields
    path "versions/*"
    path "params.json"
    val wf_version
    
    output:
    path "wf-bam-merge-*.html"
    
    script:
    String report_name = analysis_group ? \
        "wf-bam-merge-$analysis_group-report.html" : "wf-bam-merge-report.html"
    String metadata_json = new JsonBuilder(metadata).toPrettyString().replaceAll("'", "'\\\\''")
    String group_arg = analysis_group ? "--analysis_group $analysis_group" : ""
    String merged_args = merged_bams ? "--merged_bams $merged_bams" : ""
    String client_fields_args = client_fields.name == OPTIONAL_FILE.name ? "" : "--client_fields $client_fields"
    
    """
    echo '${metadata_json}' > metadata.json
    workflow-glue report $report_name \\
        $group_arg \\
        --versions versions \\
        $merged_args \\
        $client_fields_args \\
        --params params.json \\
        --metadata metadata.json \\
        --wf_version $wf_version
    """
}

process publish {
    label "wf_bam_merge"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    
    input:
    tuple path(fname), val(dirname)
    
    output:
    path fname
    
    """
    echo "Publishing output files"
    """
}
process renameBamFiles {
    label "wf_bam_merge"
    
    input:
    tuple val(meta), path(bam), path(bai), path(stats)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), path(stats), emit: renamed_bam
    
    script:
    // Build custom filename based on parameters
    def sample_name = params.sample ?: meta.alias ?: meta.id ?: "sample"
    def custom_prefix = params.output_prefix ?: sample_name
    def alignment_suffix = params.include_alignment_status ? (params.alignment_status == 'aligned' ? ".aligned" : ".unaligned") : ""
    def date_suffix = params.include_date ? ".${new Date().format('yyyyMMdd')}" : ""
    def output_name = "${custom_prefix}${alignment_suffix}${date_suffix}.merged.bam"
    
    """
    echo "Renaming BAM file from ${bam} to ${output_name}"
    echo "Original meta.alias: ${meta.alias}"
    echo "Using sample name: ${sample_name}"
    
    # Copy/rename the BAM file with custom name
    cp ${bam} ${output_name}
    
    # Copy/rename the BAI file with custom name
    if [ -f "${bai}" ]; then
        cp ${bai} ${output_name}.bai
    else
        echo "No BAI file found, creating index..."
        samtools index ${output_name}
    fi
    
    echo "Successfully renamed to: ${output_name}"
    ls -la *.bam*
    """
}

process collectIngressResultsInDir {
    label "wf_bam_merge"
    
    input:
    tuple val(meta),
          path(reads, stageAs: "reads/*"),
          path(index, stageAs: "index/*"),
          path(stats, stageAs: "stats/*")
    
    output:
    path "out/*"
    
    script:
    String outdir = "out/${meta["alias"].replaceAll("'", "'\\\\''")}"
    String metaJson = new JsonBuilder(meta).toPrettyString().replaceAll("'", "'\\\\''")
    String reads = reads.fileName.name == OPTIONAL_FILE.name ? "" : reads
    String index = index.fileName.name == OPTIONAL_FILE.name ? "" : index
    String stats = stats.fileName.name == OPTIONAL_FILE.name ? "" : stats
    
    """
    mkdir -p '$outdir'
    echo '$metaJson' > metamap.json
    mv metamap.json $reads $stats $index '$outdir'
    """
}

// workflow module
workflow pipeline {
    take:
    bam_data
    
    main:
    client_fields = params.client_fields && file(params.client_fields).exists() ? file(params.client_fields) : OPTIONAL_FILE
        
    
    software_versions = getVersions()
    workflow_params = getParams()
    
    // CRITICAL FIX: Override sample names with --sample parameter and rename BAM files
    fixed_bam_data = bam_data
        .map { meta, bam, bai, stats ->
            def new_meta = [:]
            new_meta.putAll(meta) // Copy all existing metadata
            
            // Override with --sample parameter if provided
            if (params.sample) {
                new_meta.alias = params.sample
                new_meta.id = params.sample
                new_meta.barcode = params.sample
            }
            
            // Debug output
            log.info "POST-INGRESS: Original alias='${meta.alias}', New alias='${new_meta.alias}', BAM staged as='${bam.name}'"
            
            [new_meta, bam, bai, stats]
        }
    
    // Rename BAM files with custom naming scheme
    renamed_bam_data = renameBamFiles(fixed_bam_data)

     // Group renamed BAM files by sample for merging
    grouped_bams = renamed_bam_data.renamed_bam
        .map { meta, bam, bai, stats -> 
            // Group by alias, collect BAMs and indices separately
            [meta.alias, meta, bam, bai, stats]
        }
        .groupTuple(by: 0)
           .map { alias, metas, bams, bais, stats ->
            // Take the first meta (they should be the same for grouped samples)
            // Filter out null values and create proper collections
            def clean_bams = bams.findAll { it != null }
            def clean_bais = bais.findAll { it != null }
            [metas[0], clean_bams, clean_bais, stats[0]]
        }
    
    // Split into single BAMs and multiple BAMs for different handling
    grouped_bams_branched = grouped_bams.branch { meta, bams, bais, stats ->
        single_bam: bams.size() == 1
        multiple_bams: bams.size() > 1
    }
    
    // For single BAMs, just publish them directly (they're already renamed)
    single_bam_results = grouped_bams_branched.single_bam
        .map { meta, bams, bais, stats ->
            [meta, bams[0], bais[0]]
        }
    
    // Only merge if multiple BAMs
    merge_candidates = grouped_bams_branched.multiple_bams
    
     // Sort BAMs if needed
    if (params.sort_before_merge) {
        // Flatten for individual sorting, then regroup
        sorted_bams = merge_candidates
            .flatMap { meta, bams, bais, stats ->
                bams.withIndex().collect { bam, idx -> 
                    [meta + [sort_index: idx], bam, bais[idx] ?: OPTIONAL_FILE]
                }
            }
            | SAMTOOLS_SORT
            | map { meta, sorted_bam, sorted_bai ->
                // Remove sort_index and regroup by alias
                def clean_meta = meta.findAll { it.key != 'sort_index' }
                [clean_meta.alias, clean_meta, sorted_bam, sorted_bai]
            }
            | groupTuple(by: 0)
            | map { alias, metas, sorted_bams, sorted_bais ->
                [metas[0], sorted_bams, sorted_bais, null]
            }
        
        merge_input = sorted_bams
    } else {
        merge_input = merge_candidates
    }
    
      // Merge BAM files using appropriate method based on alignment status
    if (params.alignment_status == 'aligned') {
        merged_results = SAMTOOLS_MERGE_ALIGNED(merge_input)
    } else {
        merged_results = SAMTOOLS_CAT_UNALIGNED(merge_input)
    }
    
    // Combine single BAMs and merged BAMs for final results
    all_final_bams = single_bam_results.mix(merged_results.merged_bam)
    
    // Prepare data for reporting
    for_report = merged_results.merged_bam
        | map { meta, merged_bam, merged_bai ->
            [meta.analysis_group ?: null, meta, merged_bam]
        }
        | groupTuple(by: 0)
        | map { analysis_group, metas, merged_bams ->
            [analysis_group, metas, merged_bams]
        }
    
    // Generate report
    report = makeReport(
        for_report,
        client_fields,
        software_versions,
        workflow_params,
        workflow.manifest.version
    )
    
    // Collect ingress results (following EPI2ME pattern) - use renamed BAM data
    renamed_bam_data.renamed_bam
        | map { meta, bam, bai, stats ->
            [meta, bam ?: OPTIONAL_FILE, bai ?: OPTIONAL_FILE, stats ?: OPTIONAL_FILE]
        }
        | collectIngressResultsInDir

    emit:
    final_bams = all_final_bams
    ingress_results = collectIngressResultsInDir.out
    report
    telemetry = workflow_params
}

// entrypoint workflow (following EPI2ME pattern)
WorkflowMain.initialise(workflow, params, log)

workflow {
    Pinguscript.ping_start(nextflow, workflow, params)
    
    // Use EPI2ME xam_ingress for BAM input handling
    def samples = xam_ingress([
        "input": params.bam,
        "sample": params.sample,
        "sample_sheet": params.sample_sheet,
        "keep_unaligned": params.wf.keep_unaligned,
        "stats": params.wf.bamstats,
        "return_fastq": false,  // We want BAM output
        "fastq_chunk": null,
        "per_read_stats": params.wf.per_read_stats,
        "allow_multiple_basecall_models": params.wf.allow_multiple_basecall_models,
    ])
    
    // Run the pipeline
    pipeline(samples)
    
    // Publish results (following EPI2ME pattern)
    ch_to_publish = pipeline.out.ingress_results
        | map { [it, params.ingress_results_dir] }
    
    ch_to_publish | publish
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}

workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}