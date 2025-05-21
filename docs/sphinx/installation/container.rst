Container-Based Installation
============================

This guide covers installing and running METADA using containers, which provides a consistent environment across different platforms.

Singularity/Apptainer Installation
----------------------------------

Prerequisites
~~~~~~~~~~~~~

1. Install Singularity/Apptainer on your system:

   Linux:
   
   .. code-block:: bash

      # For Ubuntu/Debian
      sudo apt update
      sudo apt install singularity-container

      # For RHEL/CentOS
      sudo yum install singularity

   macOS:
   
   .. code-block:: bash

      # Using Homebrew
      brew install singularity

   Windows:
   
   Use Windows Subsystem for Linux (WSL2) and follow Linux instructions.

Using Pre-built Container
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Pull the METADA container:

   .. code-block:: bash

      singularity pull metada.sif oras://ghcr.io/metada/metada:latest

2. Run the container:

   .. code-block:: bash

      # Linux/WSL with GPU support
      singularity shell --nv metada.sif

      # macOS or without GPU
      singularity shell metada.sif

Building Custom Container
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create a definition file ``metada.def``:

   .. literalinclude:: ../../../environments/container/singularity/metada.def
      :language: singularity
      :caption: metada.def

2. Build the container:

   .. code-block:: bash

      sudo singularity build metada.sif metada.def

Docker Installation
-------------------

Prerequisites
~~~~~~~~~~~~~

1. Install Docker:
   
   - `Docker Desktop for Windows <https://docs.docker.com/desktop/install/windows-install/>`_
   - `Docker Desktop for macOS <https://docs.docker.com/desktop/install/mac-install/>`_
   - `Docker Engine for Linux <https://docs.docker.com/engine/install/>`_

2. For GPU support:
   
   - Linux: Install `NVIDIA Container Toolkit <https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html>`_
   - Windows: Enable WSL2 backend in Docker Desktop

Using Pre-built Container
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Pull the METADA container:

   .. code-block:: bash

      docker pull ghcr.io/metada/metada:latest

2. Run the container:

   .. code-block:: bash

      # With GPU support
      docker run --gpus all -it ghcr.io/metada/metada:latest

      # Without GPU
      docker run -it ghcr.io/metada/metada:latest

Building Custom Container
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create a Dockerfile:

   .. literalinclude:: ../../../environments/container/docker/Dockerfile
      :language: dockerfile
      :caption: Dockerfile


2. Build the container:

   .. code-block:: bash

      docker build -t metada:custom .

Development Workflow
--------------------

VS Code Integration
~~~~~~~~~~~~~~~~~~~

1. Install required VS Code extensions:
   
   - Remote Development
   - Docker (for Docker workflow)

2. Configure container development:

   For Docker:
   
   .. code-block:: bash

      # Create a .devcontainer directory
      mkdir -p .devcontainer
      
      # Create a basic devcontainer configuration
      cat > .devcontainer/devcontainer.json << 'EOF'
      {
          "name": "METADA Development",
          "image": "ghcr.io/metada/metada:latest",
          "runArgs": ["--gpus", "all"],
          "customizations": {
              "vscode": {
                  "extensions": [
                      "ms-vscode.cpptools",
                      "ms-vscode.cmake-tools",
                      "twxs.cmake"
                  ]
              }
          },
          "remoteUser": "vscode"
      }
      EOF

3. Mount source code and build:

   Singularity:
   
   .. code-block:: bash

      singularity shell --nv -B /path/to/source:/opt/metada metada.sif

   Docker:
   
   .. code-block:: bash

      docker run --gpus all -v /path/to/source:/opt/metada -it metada:latest

Building and Testing
~~~~~~~~~~~~~~~~~~~~

Inside the container:

.. code-block:: bash

   cd /opt/metada
   cmake -S . -B build
   cmake --build build -j$(nproc)
   ctest --test-dir build --output-on-failure

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

- **GPU not available**:
  - Verify NVIDIA drivers are installed
  - Check container GPU support flags (``--nv`` for Singularity, ``--gpus all`` for Docker)
  - Ensure NVIDIA Container Toolkit is installed (Docker on Linux)

- **Permission denied**:
  - Check bind mount permissions
  - Verify user permissions in container
  - Use appropriate flags for volume mounts

- **Build failures**:
  - Verify all dependencies are included in container
  - Check compiler compatibility
  - Ensure sufficient resources (memory/disk space)

- **Performance issues**:
  - Monitor GPU utilization
  - Check memory allocation
  - Verify container resource limits 