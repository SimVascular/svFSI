# Docker-svFSI
This Dockerfile will build [`svFSI`](https://github.com/SimVascular/svFSI) executable from the most recent source code in the main repository. This procedure has been successfully tested on MacOS Big Sur, Ubuntu 18.04 and Windows 10 with WSL 2. Assuming you already have [Docker](https://docs.docker.com/get-docker/) installed, please follow the steps below to run `svFSI`.

1. Build Docker image. In the current directory (Path_to_svFSI/Docker), run the following command.

   ```bash
   docker build -t svfsi-image .
   ```

   This may take a while. Afterwards, run the command `docker images`, and you should see `svfsi-image`.

3. Download the examples.

   ```bash
   git clone https://github.com/SimVascular/svFSI-Tests
   ```

4. Run the container in interactive mode.

   ```bash
   docker container run --cap-add=SYS_PTRACE -v "$PWD"/svFSI-Tests:/home/test/svFSI-Tests -it --rm --name svfsi-demo svfsi-image
   ```

   This will open a shell prompt and you can proceed as usual. Here, `--cap-add=SYS_PTRACE` fixes a known [issue](https://github.com/open-mpi/ompi/issues/4948) of running openmpi in Docker.

5. Let's take `04-fluid/06-channel-flow-2D` for example. In the shell prompt, run the following commands to generate the simulation results.

   ```bash
   cd svFSI-Tests/04-fluid/06-channel-flow-2D && \
   mpiexec -n 4 svFSI ./svFSI_Taylor-Hood.inp
   ```

   The results will be stored in `4-procs` in vtu format, and can be viewed with [Paraview](https://www.paraview.org).

6. After all tests are done, execute the following commands to exit the docker and delete the image.

   ```bash
   exit && \
   docker rmi <IMAGE ID>
   ```



### Known Issues

`svFSI` built with this Dockerfile won't work with any example that requires Trilinos. Trilinos takes too long to build within the Docker image, and we encourage any user that needs it to build `svFSI` from source. Please report any other issue through the GitHub page.