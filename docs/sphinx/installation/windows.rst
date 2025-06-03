Windows Installation Guide
==========================

Setting up Development Environment with MSYS2
---------------------------------------------

Prerequisites
~~~~~~~~~~~~~

1. **NVIDIA GPU** with CUDA support (optional, for GPU acceleration)
2. **MSYS2** (required)
   - Download from `MSYS2 website <https://www.msys2.org/>`_
   - Run the installer and follow the installation wizard

Installing Required Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Open "MSYS2 MinGW64" terminal and update the system:

   .. code-block:: bash

      pacman -Syu

2. Install required development packages:

   .. code-block:: bash

      # Core development tools
      pacman -S mingw-w64-x86_64-gcc
      pacman -S mingw-w64-x86_64-gcc-fortran
      pacman -S mingw-w64-x86_64-cmake
      pacman -S mingw-w64-x86_64-ninja

      # Python and dependencies
      pacman -S mingw-w64-x86_64-python
      pacman -S mingw-w64-x86_64-python-pip
      pacman -S mingw-w64-x86_64-python-numpy

      # Testing and development tools
      pacman -S mingw-w64-x86_64-gtest
      pacman -S mingw-w64-x86_64-glog
      pacman -S mingw-w64-x86_64-clang-tools-extra
      pacman -S mingw-w64-x86_64-lcov

      # Configuration backends
      pacman -S mingw-w64-x86_64-yaml-cpp
      pacman -S mingw-w64-x86_64-nlohmann-json

      # Logging library (install from source)
      git clone --depth 1 https://github.com/ng-log/ng-log.git
      cd ng-log
      mkdir build
      cd build
      cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/mingw64
      ninja
      ninja install
      cd ../..

IDE Setup
~~~~~~~~~

1. Install VS Code from the `official website <https://code.visualstudio.com/>`_
2. Install required extensions:

   .. code-block:: bash

      code --install-extension ms-vscode.cpptools
      code --install-extension ms-vscode.cmake-tools
      code --install-extension twxs.cmake

Building METADA
~~~~~~~~~~~~~~~

1. Configure the project:

   .. code-block:: bash

      mkdir build
      cd build
      cmake -G Ninja \
        -DPython3_ROOT_DIR=/mingw64 \
        -DPython3_EXECUTABLE=/mingw64/bin/python3 \
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

Note: Python tests are disabled on Windows with MSYS2.

Troubleshooting
~~~~~~~~~~~~~~~

Common issues and solutions:

- **CMake not found**: Ensure MSYS2's MinGW64 bin directory (C:/msys64/mingw64/bin) is in your system PATH
- **Build errors**: Run ``pacman -Syu`` to ensure all packages are up to date
- **CUDA errors**: Install CUDA Toolkit from NVIDIA's website if GPU support is needed 