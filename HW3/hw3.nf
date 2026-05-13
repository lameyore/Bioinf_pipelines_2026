#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sra_id    = null
params.reads     = null
params.reference = null
params.outdir    = "results"

params.min_base_quality  = 20
params.min_map_quality   = 20
params.min_depth         = 5
params.ploidy            = 2


process FETCH_SRA {
    tag "$sra_id"
    conda 'bioconda::sra-tools=3.0.10'
    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}_{1,2}.fastq.gz")

    script:
    """
    prefetch $sra_id
    fasterq-dump --threads ${task.cpus} --split-files $sra_id
    pigz -p ${task.cpus} ${sra_id}_1.fastq ${sra_id}_2.fastq 2>/dev/null || \
        gzip ${sra_id}_1.fastq ${sra_id}_2.fastq
    """
}

process FASTQC_RAW {
    tag "$sample_id"
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    def files = reads instanceof List ? reads.join(' ') : reads
    """
    fastqc -t ${task.cpus} -q ${files}
    """
}

process TRIMMOMATIC {
    tag "$sample_id"
    conda 'bioconda::trimmomatic=0.39'
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: "*_p.fastq.gz"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_p.fastq.gz"), emit: reads
    path "${sample_id}_trimmomatic.log",                           emit: log

    script:
    """
    trimmomatic PE -threads ${task.cpus} \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_p.fastq.gz ${sample_id}_R1_u.fastq.gz \\
        ${sample_id}_R2_p.fastq.gz ${sample_id}_R2_u.fastq.gz \\
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \\
        2> ${sample_id}_trimmomatic.log
    """
}

process FASTQC_TRIMMED {
    tag "$sample_id"
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    def files = reads instanceof List ? reads.join(' ') : reads
    """
    fastqc -t ${task.cpus} -q ${files}
    """
}

process SPADES {
    tag "$sample_id"
    conda 'bioconda::spades=3.15.5 conda-forge::python=3.10'
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    spades.py --careful -t ${task.cpus} \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -o spades_out
    cp spades_out/contigs.fasta ${sample_id}_contigs.fasta
    """
}

process MAP_READS {
    tag "$sample_id"
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.19'
    publishDir "${params.outdir}/mapped", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    bwa index $reference

    bwa mem -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        $reference ${reads[0]} ${reads[1]} \\
    | samtools view -@ ${task.cpus} -bS -F 4 - \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}

process PLOT_COVERAGE {
    tag "$sample_id"
    conda 'bioconda::samtools=1.19 conda-forge::python=3.10 conda-forge::matplotlib=3.8.0'
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_coverage.png"), path("${sample_id}_depth.txt")

    script:
    """
    samtools depth -a $bam > ${sample_id}_depth.txt

    cat <<'EOF' > plot_coverage.py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

positions, depths = [], []
with open("${sample_id}_depth.txt") as f:
    for line in f:
        parts = line.strip().split()
        positions.append(int(parts[1]))
        depths.append(int(parts[2]))

depths = np.array(depths)

window = max(1, len(depths) // 500)
rolling = np.convolve(depths, np.ones(window)/window, mode='same')

fig, ax = plt.subplots(figsize=(14, 4))
ax.fill_between(positions, depths, alpha=0.25, color='#4C72B0', label='Per-base depth')
ax.plot(positions, rolling, color='#C44E52', linewidth=1.2,
        label=f'{window}-bp rolling mean')
ax.axhline(10, color='green',  linestyle='--', linewidth=0.8, alpha=0.7, label='10x')
ax.axhline(30, color='orange', linestyle=':',  linewidth=0.8, alpha=0.7, label='30x')

mean_d   = depths.mean()
pct_cov  = (depths > 0).sum() / len(depths) * 100
ax.set_title(f'Coverage – ${sample_id}   '
             f'(mean: {mean_d:.1f}x, covered: {pct_cov:.1f}%)', fontsize=11)
ax.set_xlabel('Position (bp)')
ax.set_ylabel('Depth')
ax.legend(fontsize=8, loc='upper right')
ax.spines[['top','right']].set_visible(False)
plt.tight_layout()
plt.savefig('${sample_id}_coverage.png', dpi=200)
print(f"Saved ${sample_id}_coverage.png  (mean depth: {mean_d:.1f}x, {pct_cov:.1f}% covered)")
EOF

    python3 plot_coverage.py
    """
}



process BCFTOOLS_CALL {
    tag "$sample_id"
    conda 'bioconda::bcftools=1.19 bioconda::samtools=1.19'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.vcf.gz.tbi"), emit: tbi
    path "${sample_id}_bcftools_stats.txt",                 emit: stats

    script:
    """
    # Index the reference if .fai is not present
    [ -f ${reference}.fai ] || samtools faidx $reference

    # mpileup → call pipeline (equivalent to nf-core BCFTOOLS_MPILEUP + BCFTOOLS_CALL)
    bcftools mpileup \\
        --fasta-ref $reference \\
        --min-BQ ${params.min_base_quality} \\
        --min-MQ ${params.min_map_quality} \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD \\
        --output-type u \\
        $bam \\
    | bcftools call \\
        --multiallelic-caller \\
        --ploidy ${params.ploidy} \\
        --variants-only \\
        --output-type z \\
        --output ${sample_id}.vcf.gz

    bcftools index --tbi ${sample_id}.vcf.gz

    # Basic stats (mirrors nf-core BCFTOOLS_STATS)
    bcftools stats ${sample_id}.vcf.gz > ${sample_id}_bcftools_stats.txt
    """
}

process FILTER_VARIANTS {
    tag "$sample_id"
    conda 'bioconda::bcftools=1.19'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz.tbi"), emit: tbi

    script:
    """
    bcftools filter \\
        --include 'INFO/DP >= ${params.min_depth} && QUAL >= 20' \\
        --output-type z \\
        --output ${sample_id}.filtered.vcf.gz \\
        $vcf

    bcftools index --tbi ${sample_id}.filtered.vcf.gz
    """
}


workflow QC_TRIM {
    take:
    raw_reads_ch   // tuple val(sample_id), path(reads)

    main:
    FASTQC_RAW(raw_reads_ch)
    trimmed_ch = TRIMMOMATIC(raw_reads_ch).reads
    FASTQC_TRIMMED(trimmed_ch)

    emit:
    reads = trimmed_ch
}

workflow ASSEMBLE_OR_MAP {
    take:
    trimmed_ch     
    ref_supplied  

    main:
    if (ref_supplied) {
        ref_ch        = Channel.fromPath(params.reference, checkIfExists: true).first()
        mapping_in_ch = trimmed_ch.combine(ref_ch)
    } else {
        assembly_ch   = SPADES(trimmed_ch)
        mapping_in_ch = trimmed_ch.join(assembly_ch)
    }

    bam_ch = MAP_READS(mapping_in_ch)
    PLOT_COVERAGE(bam_ch)

    emit:
    bam = bam_ch
}

workflow VARIANT_CALLING {
    take:
    bam_ch   // tuple val(sample_id), path(bam), path(bai)

    main:
    ref_ch       = Channel.fromPath(params.reference, checkIfExists: true).first()

    // Attach reference to every BAM tuple
    bam_ref_ch   = bam_ch.combine(ref_ch)

    raw_vcf_ch   = BCFTOOLS_CALL(bam_ref_ch)
    FILTER_VARIANTS(raw_vcf_ch.vcf.join(raw_vcf_ch.tbi))

    emit:
    vcf = FILTER_VARIANTS.out.vcf
}



workflow {


    if (params.sra_id) {
        raw_reads_ch = FETCH_SRA(Channel.of(params.sra_id))
    } else if (params.reads) {
        raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        error "Provide --sra_id <accession>  or  --reads \"data/*_{1,2}.fastq.gz\""
    }

    trimmed_ch = QC_TRIM(raw_reads_ch).reads


    bam_ch = ASSEMBLE_OR_MAP(trimmed_ch, params.reference as Boolean).bam


    if (!params.reference) {
        log.warn "No --reference provided: variant calling will be skipped (de-novo assembly used as reference for mapping only)."
    } else {
        VARIANT_CALLING(bam_ch)
    }
}
