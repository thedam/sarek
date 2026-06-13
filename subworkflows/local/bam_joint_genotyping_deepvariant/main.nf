//
// JOINT GENOTYPING OF DEEPVARIANT gVCFs
//
// Merge the per-sample DeepVariant gVCFs of the whole cohort and perform
// joint genotyping with GLnexus, then convert the resulting multi-sample BCF
// to a bgzipped, tabix-indexed VCF ready for annotation.
//
// GLnexus is the joint genotyper recommended for DeepVariant gVCFs (it ships a
// dedicated `DeepVariant*` configuration preset). The preset is supplied via
// `ext.args` in conf/modules/joint_genotype.config.
//

include { GLNEXUS                              } from '../../../modules/nf-core/glnexus/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_JOINT } from '../../../modules/nf-core/bcftools/view/main'

workflow BAM_JOINT_GENOTYPING_DEEPVARIANT {
    take:
    gvcf // channel: [mandatory] [ meta, gvcf ] one per cohort sample
    bed  // channel: [optional]  [ meta, bed ] regions to restrict joint genotyping to, [ [id:'no_region'], [] ] for whole genome

    main:
    versions = Channel.empty()

    // Gather every sample gVCF into a single list and assign a cohort-level meta.
    // GLnexus reads all gVCFs in one pass and emits a single multi-sample BCF.
    glnexus_input = gvcf
        .map { _meta, gvcf -> gvcf }
        .collect()
        .map { gvcfs -> [ [ id:'joint_genotyping' ], gvcfs, [] ] }

    GLNEXUS(glnexus_input, bed)

    // Convert the multi-sample BCF to bgzipped + tabix-indexed VCF for downstream annotation.
    bcftools_input = GLNEXUS.out.bcf.map { meta, bcf -> [ meta, bcf, [] ] }

    BCFTOOLS_VIEW_JOINT(bcftools_input, [], [], [])

    // Rework meta for variantscalled.csv, post-variantcalling and annotation tools.
    genotype_vcf = BCFTOOLS_VIEW_JOINT.out.vcf
        .map { _meta, vcf -> [ [ id:'joint_genotyping', patient:'all_samples', variantcaller:'deepvariant' ], vcf ] }

    genotype_index = BCFTOOLS_VIEW_JOINT.out.tbi
        .map { _meta, tbi -> [ [ id:'joint_genotyping', patient:'all_samples', variantcaller:'deepvariant' ], tbi ] }

    // GLnexus reports its version through the `versions` channel topic, collected centrally.
    versions = versions.mix(BCFTOOLS_VIEW_JOINT.out.versions)

    emit:
    genotype_vcf   // channel: [ val(meta), [ vcf ] ]
    genotype_index // channel: [ val(meta), [ tbi ] ]

    versions       // channel: [ versions.yml ]
}
