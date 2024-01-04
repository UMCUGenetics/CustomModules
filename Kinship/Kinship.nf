process Kinship {
    tag {"Kinship ${analysis_id}"}
    label 'Kinship'
    container = 'docker.io/umcugenbioinf/kinship:1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(analysis_id), path(vcf_file), path(vcf_index))
        path(ped_file)

    output:
        tuple(val(analysis_id), path("${analysis_id}.kinship"), path("${analysis_id}.kinship_check.out"))

    script:
        """
        vcftools --vcf ${vcf_file} --plink
        plink --file out --make-bed --noweb
        king -b plink.bed --kinship
        cp king.kin0 ${analysis_id}.kinship
        python ${baseDir}/CustomModules/Kinship/check_kinship.py ${analysis_id}.kinship ${ped_file} > ${analysis_id}.kinship_check.out
        """
}