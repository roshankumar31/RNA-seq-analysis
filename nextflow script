// Define the workflow script
params.reads = "/home/vignesh/roshan/nextflow/rnaseq/SRR1634078.fastq" 
params.skip_fastqc = false 
params.with_umi = true 
params.skip_umi_extract = false 
params.umi_discard_read = 0 
params.skip_trimming = false 
params.adapter_fasta = "/path/to/adapter.fasta" 
params.save_trimmed_fail = false 
params.save_merged = false 
params.min_trimmed_reads = 100 

// Include modules
include { FASTQC as FASTQC_RAW  } from '/home/vignesh/roshan/nextflow/rnaseq/fastqc.git'
include { FASTQC as FASTQC_TRIM } from '/home/vignesh/roshan/nextflow/rnaseq/fastqc.git'
include { UMITOOLS_EXTRACT      } from '/home/vignesh/roshan/nextflow/rnaseq/umitools.git'
include { FASTP                 } from '/home/vignesh/roshan/nextflow/rnaseq/fastp.git'

// Define function to parse FASTP JSON output
import groovy.json.JsonSlurper

def getFastpReadsAfterFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toLong()
}

// Define workflow
workflow FASTQ_FASTQC_UMITOOLS_FASTP {
    take:
    reads
    skip_fastqc
    with_umi
    skip_umi_extract
    umi_discard_read
    skip_trimming
    adapter_fasta
    save_trimmed_fail
    save_merged
    min_trimmed_reads

    main:
    ch_versions = Channel.empty()
    fastqc_raw_html = Channel.empty()
    fastqc_raw_zip  = Channel.empty()
    if (!skip_fastqc) {
        FASTQC_RAW (
            reads
        )
        fastqc_raw_html = FASTQC_RAW.out.html
        fastqc_raw_zip  = FASTQC_RAW.out.zip
        ch_versions     = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    umi_reads = reads
    umi_log   = Channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT (
            reads
        )
        umi_reads   = UMITOOLS_EXTRACT.out.reads
        umi_log     = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1,2]) {
            UMITOOLS_EXTRACT
                .out
                .reads
                .map {
                    meta, reads ->
                        meta.single_end ? [ meta, reads ] : [ meta + [single_end: true], reads[umi_discard_read % 2] ]
                }
                .set { umi_reads }
        }
    }

    trim_reads        = umi_reads
    trim_json         = Channel.empty()
    trim_html         = Channel.empty()
    trim_log          = Channel.empty()
    trim_reads_fail   = Channel.empty()
    trim_reads_merged = Channel.empty()
    fastqc_trim_html  = Channel.empty()
    fastqc_trim_zip   = Channel.empty()
    trim_read_count   = Channel.empty()
    if (!skip_trimming) {
        FASTP (
            umi_reads,
            adapter_fasta,
            save_trimmed_fail,
            save_merged
        )
        trim_json         = FASTP.out.json
        trim_html         = FASTP.out.html
        trim_log          = FASTP.out.log
        trim_reads_fail   = FASTP.out.reads_fail
        trim_reads_merged = FASTP.out.reads_merged
        ch_versions       = ch_versions.mix(FASTP.out.versions.first())

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        FASTP
            .out
            .reads
            .join(trim_json)
            .map { meta, reads, json -> [ meta, reads, getFastpReadsAfterFiltering(json) ] }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { meta, reads, num_reads -> num_reads >= min_trimmed_reads.toLong() }
            .map { meta, reads, num_reads -> [ meta, reads ] }
            .set { trim_reads }

        ch_num_trimmed_reads
            .map { meta, reads, num_reads -> [ meta, num_reads ] }
            .set { trim_read_count }

        if (!skip_fastqc) {
            FASTQC_TRIM (
                trim_reads
            )
            fastqc_trim_html = FASTQC_TRIM.out.html
            fastqc_trim_zip  = FASTQC_TRIM.out.zip
            ch_versions      = ch_versions.mix(FASTQC_TRIM.out.versions.first())
        }
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_raw_html    // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip     // channel: [ val(meta), [ zip ] ]

    umi_log            // channel: [ val(meta), [ log ] ]

    trim_json          // channel: [ val(meta), [ json ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_log           // channel: [ val(meta), [ log ] ]
    trim_reads_fail    // channel: [ val(meta), [ fastq.gz ] ]
    trim_reads_merged  // channel: [ val(meta), [ fastq.gz ] ]
    trim_read_count    // channel: [ val(meta), val(count) ]

    fastqc_trim_html   // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip    // channel: [ val(meta), [ zip ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}

// Run the workflow
