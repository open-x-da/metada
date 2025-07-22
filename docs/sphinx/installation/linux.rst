Linux Installation Guide
========================

Setting up Development Environment using Spack
----------------------------------------------

Prerequisites
~~~~~~~~~~~~~

Before starting, ensure you have the following:

1. **NVIDIA GPU** with CUDA support (for GPU acceleration)
2. **NVIDIA Driver** installed on your system (if using GPU acceleration)
   - Check with ``nvidia-smi`` command
   - If not installed, follow your distribution's instructions for NVIDIA driver installation

Installing Spack
~~~~~~~~~~~~~~~~

1. Clone Spack repository:

   .. code-block:: bash

      git clone -c feature.manyFiles=true https://github.com/spack/spack.git
      cd spack
      git checkout releases/v0.21

2. Set up Spack environment:

   .. code-block:: bash

      # Add to your ~/.bashrc or ~/.zshrc
      export SPACK_ROOT=/path/to/spack
      source $SPACK_ROOT/share/spack/setup-env.sh

3. Verify installation:

   .. code-block:: bash

      spack --version
      spack compiler find

Setting up METADA Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create a new Spack environment:

   .. code-block:: bash

      spack env create metada
      spack env activate metada

2. Add required packages:

   .. code-block:: bash

      # Core development tools
      spack add cmake@3.30: ninja gcc
      # GPU acceleration support
      spack add cuda@12.1.0
      # Scientific computing packages
      spack add eigen boost
      # Configuration packages
      spack add yaml-cpp nlohmann-json
      # Testing and development tools
      spack add googletest
      spack add lcov

3. Install ng-log from source (since it's not available in Spack):

   .. code-block:: bash

      git clone --depth 1 https://github.com/ng-log/ng-log.git
      cd ng-log
      mkdir build && cd build
      cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SPACK_ROOT/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/ng-log-master
      ninja
      ninja install
      cd ../..

4. Configure environment:

   .. literalinclude:: ../../../environments/spack/spack.yaml
      :language: yaml
      :caption: spack.yaml

5. Install the environment:

   .. code-block:: bash

      spack install

IDE Setup
~~~~~~~~~

1. Install VS Code:
   
   .. code-block:: bash

      # For Debian/Ubuntu
      sudo apt-get install code
      # For other distributions, download from code.visualstudio.com

2. Install extensions:

   .. code-block:: bash

      code --install-extension ms-vscode.cpptools
      code --install-extension ms-vscode.cmake-tools
      code --install-extension twxs.cmake

Building METADA
~~~~~~~~~~~~~~~

1. Build the project:

   .. code-block:: bash

      cmake -S . -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_CXX_STANDARD=17 \
        -DCMAKE_PREFIX_PATH=$SPACK_ROOT/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/ng-log-master
      cmake --build build

Running Tests
~~~~~~~~~~~~~

Execute the test suite:

.. code-block:: bash

   cd build
   ctest --output-on-failure

Troubleshooting
~~~~~~~~~~~~~~~

Common issues and solutions:

- **Spack environment issues**: Verify environment activation
- **Build failures**: Check compiler compatibility
- **Missing dependencies**: Use ``spack spec metada`` to verify package resolution
- **ng-log issues**: Ensure ng-log is properly installed and CMAKE_PREFIX_PATH includes its installation directory 