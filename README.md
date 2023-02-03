# CustomModules
UMCU Genetics Custom Modules - Genome Diagnostics.

This repository contains custom nextflow processes and their dependent/linked files (scripts, dockerfiles etc.) used as submodule in our NGS workflows.

## Set up folder structure
### Utils
Utility functions are dinstinguished and can be found in Utils.

### Designated folder
When nextflow processes or their dependent files are linked to another git repository, 
the files will be placed in a designated folder. For example:
ClarityEpp folder and https://github.com/UMCUGenetics/clarity_epp

If a set of files is not per se a utility and doesnot have a separate repository, it is allowed to create a designated folder as well.

## Docker files
Build docker image for software dependencies. 
- [Install Docker Desktop](https://docs.docker.com/desktop/mac/apple-silicon/)
    ```bash
    docker build -t organization_or_username/toolname:version -f path_to_dockerfile
    docker push organization_or_username/toolname:version
    ```
- If changes are required to the dockerfile, manually update label version.