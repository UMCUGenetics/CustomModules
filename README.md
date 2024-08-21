# CustomModules
UMCU Genetics Custom Modules - Genome Diagnostics.

This repository contains custom nextflow processes and their dependent/linked files (scripts, dockerfiles etc.) used as submodule in our NGS workflows.

## Pytests and requirements (disclaimer)
Custom modules with dependent python scripts will be covered by pytests.
We would like to trigger these pytest automatically. For example by using Github Actions.
This could be challenging with all different package versions.
We will discuss a more futureproof way together to resolve several challenges with current setup.

In the meantime...
Two types of packages are distinguished:
- Packages required to run software (prod.txt)
- Packages required to run pytest / unit testing (testing.txt)

The requirements in a CustomModule designated folder will reference the `testing.txt` and add packages required to run software.
Note: Specify the package with exact version if possible.

When adding your custom module with pytest,
1. Create a requirements.txt in your designated folder. `CustomModules/Tool/version/requirements.txt`
   1. Add reference `-r` to `CustomModules/requirements/testing.txt`
   2. Add required packages to run software.
2. Add designated requirements reference to the generic `CustomModules/requirements/prod.txt`.
3. Resolve any conflicting dependencies.
  - Select highest version of packages required to run the software
  - Select version of testing packages as defined in `testing.txt`. Unless otherwise required, than upgrade accordingly.
4. Check if pytest succeeds if triggered in main folder of CustomModules.

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
