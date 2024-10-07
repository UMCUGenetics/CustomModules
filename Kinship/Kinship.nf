process Kinship {
    tag {"Kinship ${analysis_id}"}
    label 'Kinship'
    container = 'ghcr.io/umcugenetics/kinship:1.1.1'

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
        python ${projectDir}/CustomModules/Kinship/check_kinship.py ${analysis_id}.kinship ${ped_file} --output_prefix ${analysis_id} --output_path .
        """
}
