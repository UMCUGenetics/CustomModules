def extractRawfilesFromDir(dir) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    dir = dir.tokenize().collect{"$it/*.raw"}
    Channel
    .fromPath(dir, type:'file')
    .ifEmpty { error "No raw files found in ${dir}." }
    .map { rawfiles_path ->
        def file_id = rawfiles_path.getSimpleName()
        [file_id, rawfiles_path]
    }
}

