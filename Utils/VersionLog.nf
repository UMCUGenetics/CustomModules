process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    input:
        path(git_dirs)

    output:
        path('repository_version.log')

    script:
        """
        for git_dir in ${git_dirs}
        do
            echo "\${git_dir}" >> repository_version.log
            git --git-dir=\${git_dir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        done
        """
}