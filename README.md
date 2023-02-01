# CustomModules
## Docker files
Build docker image for software dependencies. 
- [Install Docker Desktop](https://docs.docker.com/desktop/mac/apple-silicon/)
    ```bash
    docker build -t organization_or_username/toolname:version -f path_to_dockerfile
    docker push organization_or_username/toolname:version
    ```
- If changes are required to the dockerfile, manually update label version.