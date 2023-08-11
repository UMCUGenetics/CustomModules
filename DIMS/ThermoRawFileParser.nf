process ConvertRawFile {
    // Custom process to convert raw file to mzML format
    tag {"DIMS ConvertRawFile ${file_id}"}
    label 'ThermoRawFileParser_1_1_11'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(file_id, path(raw_file))

    output:
        tuple(file_id, path("${file_id}.mzML"))

    script:

        """
        source /hpc/dbg_mz/tools/mono/etc/profile
        mono /hpc/dbg_mz/tools/ThermoRawFileParser_1.1.11/ThermoRawFileParser.exe -i=${raw_file} --output=./
        """
}
