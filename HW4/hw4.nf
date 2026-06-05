#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREFETCH      } from './modules/nf-core/sratools/prefetch/main'
include { FASTERQ_DUMP  } from './modules/nf-core/sratools/fasterqdump/main'
include { FASTQC as FASTQC_RAW  } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc/main'
include { TRIMMOMATIC   } from './modules/nf-core/trimmomatic/main'
include { SPADES        } from './modules/nf-core/spades/main'
include { BWA_MEM       } from './modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT  } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH } from './modules/nf-core/samtools/depth/main'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CALL    } from './modules/nf-core/bcftools/call/main'
include { BCFTOOLS_FILTER  } from './modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_STATS   } from './modules/nf-core/bcftools/stats/main'

params.sra_id            = null
params.reads             = null
params.reference         = null
params.samplesheet       = null  
params.outdir            = "results"
params.min_base_quality  = 20
params.min_map_quality   = 20
params.min_depth         = 5
params.ploidy            = 2

process SUMMARIZE_GROUP {
    tag "${group_id}"
    publishDir "${params.outdir}/summary/${group_id}", mode: 'copy'

    input:
    tuple val(group_id), path(vcf_files)

    output:
    tuple val(group_id), path("${group_id}_summary.tsv"), emit: summary
    """
    echo -e "group\tsample\tnum_variants\tpassed_variants" > ${group_id}_summary.tsv
    for vcf in ${vcf_files}; do
        sample=\$(basename "\$vcf" .vcf.gz)
        echo -e "${group_id}\t\${sample}\t0\t0" >> ${group_id}_summary.tsv
    done
    """
    script:
    """
    echo -e "group\tsample\tnum_variants\tpassed_variants" > ${group_id}_summary.tsv

    for vcf in ${vcf_files}; do
        sample=\$(basename "\$vcf" .vcf.gz)
        total=\$(bcftools view -H "\$vcf" | wc -l)
        passed=\$(bcftools view -f PASS -H "\$vcf" | wc -l)
        echo -e "${group_id}\t\${sample}\t\${total}\t\${passed}" >> ${group_id}_summary.tsv
    done
    """
}

def parseSamplesheet(csv_path) {
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id        : row.sample_id,
                group     : row.group,
                reference : row.reference
            ]
            def reads = [ file(row.read1), file(row.read2) ]
            return tuple(meta, reads)
        }
}

workflow FETCH_SRA {
    take:
    sra_ch   

    main:
    prefetch_ch = PREFETCH(sra_ch)
    fastq_ch    = FASTERQ_DUMP(prefetch_ch)

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
     
    
        map_input = reads_ch.map { meta, reads ->
            def ref = meta.reference ? file(meta.reference) : file(params.reference)
            tuple(meta, reads, ref)
        }
        bam = BWA_MEM(map_input)
    } else {
        contigs   = SPADES(reads_ch)
        map_input = reads_ch.join(contigs).map { meta, reads, ctg -> tuple(meta, reads, ctg) }
        bam       = BWA_MEM(map_input)
    }

    sorted  = SAMTOOLS_SORT(bam)
    indexed = SAMTOOLS_INDEX(sorted)
    SAMTOOLS_DEPTH(sorted)

    emit:
    bam = sorted  
}

workflow VARIANT_CALLING {
    take:
    bam_ch   

    main:
    input_ch = bam_ch.map { meta, bam ->
        def ref = meta.reference ? file(meta.reference) : file(params.reference)
        tuple(meta, bam, ref)
    }

    mpileup  = BCFTOOLS_MPILEUP(input_ch)
    calls    = BCFTOOLS_CALL(mpileup)
    BCFTOOLS_STATS(calls)
    filtered = BCFTOOLS_FILTER(calls)

    emit:
    vcf = filtered   
}

workflow {


    if (params.samplesheet) {
        raw_reads_ch = parseSamplesheet(params.samplesheet)

    } else if (params.sra_id) {
        def meta     = [ id: params.sra_id, group: 'default', reference: params.reference ]
        raw_reads_ch = FETCH_SRA(Channel.of(tuple(meta, params.sra_id))).reads

    } else if (params.reads) {
        raw_reads_ch = Channel
            .fromFilePairs(params.reads)
            .map { sid, files ->
                def meta = [ id: sid, group: 'default', reference: params.reference ]
                tuple(meta, files)
            }

    } else {
        error "Provide --samplesheet, --sra_id, or --reads"
    }

    trimmed_ch = QC_TRIM(raw_reads_ch).reads

    branched = trimmed_ch.branch { meta, reads ->
        virus_A  : meta.group == 'virus_A'
        virus_B  : meta.group == 'virus_B'
        other    : true            // catch-all for any other group name
    }
    def ref_supplied = params.reference != null || params.samplesheet != null

    bam_A     = ASSEMBLE_OR_MAP(branched.virus_A, ref_supplied).bam
    bam_B     = ASSEMBLE_OR_MAP(branched.virus_B, ref_supplied).bam
    bam_other = ASSEMBLE_OR_MAP(branched.other,   ref_supplied).bam

    vcf_A     = VARIANT_CALLING(bam_A).vcf
    vcf_B     = VARIANT_CALLING(bam_B).vcf
    vcf_other = VARIANT_CALLING(bam_other).vcf

    grouped_vcf_ch = vcf_A
        .mix(vcf_B, vcf_other)
        .map   { meta, vcf -> tuple(meta.group, vcf) }
        .groupTuple()
    SUMMARIZE_GROUP(grouped_vcf_ch)
}
