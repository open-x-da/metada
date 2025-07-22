macOS Installation Guide
========================

.. note::
   While GPU acceleration via CUDA is not supported on macOS, METADA can still be used with CPU computations 
   or alternative GPU computing solutions.

Prerequisites
~~~~~~~~~~~~~

1. **Homebrew** package manager
2. **Development tools** and libraries

Installing Dependencies
~~~~~~~~~~~~~~~~~~~~~~~

1. Install GFortran and core dependencies:

   .. code-block:: bash

      brew update
      # Install gfortran from gcc
      brew install gcc@13
      
      # Install other dependencies
      brew install \
        cmake@3.30: \
        ninja \
        python \
        numpy \
        googletest \
        llvm \
        lcov \
        yaml-cpp \
        nlohmann-json

2. Install ng-log from source:

   .. code-block:: bash

      git clone --depth 1 https://github.com/ng-log/ng-log.git
      cd ng-log
      mkdir build && cd build
      cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
      ninja
      ninja install
      cd ../..

Building METADA
~~~~~~~~~~~~~~~

1. Configure the project:

   .. code-block:: bash

      mkdir build
      cd build
      GFORTRAN_PATH=$(brew --prefix gcc@13)/bin/gfortran-13
      cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_CXX_STANDARD=17 \
        -DCMAKE_Fortran_COMPILER=${GFORTRAN_PATH} \
        ..

2. Build the project:

   .. code-block:: bash

      cmake --build .

Running Tests
~~~~~~~~~~~~~

Execute the test suite:

.. code-block:: bash

   cd build
   ctest --output-on-failure

Alternative 1: Using Remote Development
---------------------------------------

This approach involves developing on macOS but building and running on a remote Linux machine with GPU support.

Prerequisites
~~~~~~~~~~~~~

1. **Remote Linux machine** with NVIDIA GPU (for GPU acceleration)
2. **VS Code** on your macOS system
3. **SSH access** to the remote machine

Setup Instructions
~~~~~~~~~~~~~~~~~~

1. Install VS Code on macOS
2. Install the "Remote - SSH" extension
3. Configure SSH connection to your remote machine
4. Set up the remote environment:
   
   a. Install required packages on remote machine:
      
      .. code-block:: bash
         
         sudo apt update
         sudo apt install build-essential cmake ninja-build python3 python3-pip

   b. Install CUDA Toolkit on remote machine (see Linux Installation Guide)
   
   c. Configure VS Code Remote SSH:
      - Press ``Cmd+Shift+P``
      - Select "Remote-SSH: Connect to Host"
      - Enter your SSH connection details

5. Clone and build on remote machine:
   
   .. code-block:: bash
      
      git clone https://github.com/your-org/metada.git
      cd metada
      cmake -S . -B build
      cmake --build build -j$(nproc)

Alternative 2: Using Docker Containers
--------------------------------------

This approach uses NVIDIA Docker containers for CUDA development.

Prerequisites
~~~~~~~~~~~~~

1. **Docker Desktop** for macOS
2. **Remote machine** with NVIDIA GPU (for running containers)

Setup Instructions
~~~~~~~~~~~~~~~~~~

1. Install Docker Desktop for macOS
2. Set up remote Docker context:
   
   .. code-block:: bash
      
      # Create context for remote machine
      docker context create remote --docker "host=ssh://user@remote-host"
      # Switch to remote context
      docker context use remote

3. Configure remote machine:
   
   a. Install NVIDIA Container Toolkit on remote host
   b. Configure Docker daemon for NVIDIA runtime
   c. Verify GPU access:
      
      .. code-block:: bash
         
         docker run --gpus all nvidia/cuda:12.1.0-base nvidia-smi

4. Development workflow:
   
   a. Use VS Code with Remote-Containers extension
   b. Open project in container using provided devcontainer configuration
   c. Build and run as specified in container documentation

Alternative 3: Using Cloud Services
-----------------------------------

This approach leverages cloud GPU instances for development.

Available Options
~~~~~~~~~~~~~~~~~

1. **Google Colab**
2. **AWS SageMaker**
3. **Azure ML Studio**

Setup Instructions
~~~~~~~~~~~~~~~~~~

1. Google Colab
   
   a. Upload project notebooks to Google Drive
   b. Configure GPU runtime:
      - Runtime → Change runtime type → GPU
   c. Install required packages:
      
      .. code-block:: bash
         
         !pip install cmake ninja
         !git clone https://github.com/your-org/metada.git

2. AWS SageMaker
   
   a. Launch SageMaker notebook instance with GPU
   b. Choose ML instance type with NVIDIA GPU
   c. Use provided container image or custom container
   d. Configure Git repository integration

3. Azure ML Studio
   
   a. Create compute instance with GPU
   b. Use provided Jupyter notebooks
   c. Configure development environment:
      - Install required extensions
      - Set up Git integration
      - Configure GPU compute targets 