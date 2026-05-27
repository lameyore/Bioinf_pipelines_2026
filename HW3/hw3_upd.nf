#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREFETCH } from './modules/nf-core/sratools/prefetch/main'
include { FASTERQ_DUMP } from './modules/nf-core/sratools/fasterqdump/main'
include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc/main'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main'
include { SPADES } from './modules/nf-core/spades/main'
include { BWA_MEM } from './modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH } from './modules/nf-core/samtools/depth/main'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CALL } from './modules/nf-core/bcftools/call/main'
include { BCFTOOLS_FILTER } from './modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_STATS } from './modules/nf-core/bcftools/stats/main'

params.sra_id = null
params.reads = null
params.reference = null
params.outdir = "results"
params.min_base_quality = 20
params.min_map_quality = 20
params.min_depth = 5
params.ploidy = 2

workflow FETCH_SRA {
    take:
    sra_ch

    main:
    prefetch_ch = PREFETCH(sra_ch)
    fastq_ch = FASTERQ_DUMP(prefetch_ch)

    emit:
    reads = fastq_ch
}

workflow QC_TRIM {
    take:
    raw_reads_ch

    main:
    FASTQC_RAW(raw_reads_ch)
    trimmed_ch = TRIMMOMATIC(raw_reads_ch)
    FASTQC_TRIM(trimmed_ch)

    emit:
    reads = trimmed_ch
}

workflow ASSEMBLE_OR_MAP {
    take:
    reads_ch
    ref_supplied

    main:
    if (ref_supplied) {
        ref_ch = Channel.fromPath(params.reference)
        map_input = reads_ch.combine(ref_ch)
        bam = BWA_MEM(map_input)
    } else {
        contigs = SPADES(reads_ch)
        map_input = reads_ch.combine(contigs)
        bam = BWA_MEM(map_input)
    }

    sorted = SAMTOOLS_SORT(bam)
    indexed = SAMTOOLS_INDEX(sorted)

    depth = SAMTOOLS_DEPTH(sorted)

    emit:
    bam = sorted
}

workflow VARIANT_CALLING {
    take:
    bam_ch

    main:
    ref_ch = Channel.fromPath(params.reference)

    input_ch = bam_ch.combine(ref_ch)

    mpileup = BCFTOOLS_MPILEUP(input_ch)
    calls = BCFTOOLS_CALL(mpileup)
    stats = BCFTOOLS_STATS(calls)
    filtered = BCFTOOLS_FILTER(calls)

    emit:
    vcf = filtered
}

workflow {
    if (params.sra_id) {
        raw_reads_ch = FETCH_SRA(Channel.of(params.sra_id)).reads
    } else if (params.reads) {
        raw_reads_ch = Channel.fromFilePairs(params.reads)
    } else {
        error "Provide --sra_id or --reads"
    }

    trimmed = QC_TRIM(raw_reads_ch).reads

    bam = ASSEMBLE_OR_MAP(trimmed, params.reference != null).bam

    if (params.reference) {
        VARIANT_CALLING(bam)
    }
}