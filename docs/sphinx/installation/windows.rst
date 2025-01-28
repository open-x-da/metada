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

      pacman -S mingw-w64-x86_64-gcc
      pacman -S mingw-w64-x86_64-cmake
      pacman -S mingw-w64-x86_64-ninja
      pacman -S mingw-w64-x86_64-python
      pacman -S mingw-w64-x86_64-python-numpy
      pacman -S mingw-w64-x86_64-gtest
      pacman -S mingw-w64-x86_64-glog
      pacman -S mingw-w64-x86_64-clang-tools-extra
      pacman -S mingw-w64-x86_64-lcov

3. Install optional packages:

   .. code-block:: bash

      # For CUDA support (optional)
      pacman -S mingw-w64-x86_64-cuda

      # For documentation generation (optional)
      pacman -S mingw-w64-x86_64-doxygen
      pacman -S mingw-w64-x86_64-python-sphinx
      pacman -S mingw-w64-x86_64-graphviz

Installing Visual Studio Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download VS Code from the `official website <https://code.visualstudio.com/>`_
2. Install these VS Code extensions:
   
   .. code-block:: bash

      code --install-extension ms-vscode.cpptools
      code --install-extension ms-vscode.cmake-tools
      code --install-extension twxs.cmake

Configuring VS Code
~~~~~~~~~~~~~~~~~~

1. Create or modify ``.vscode/settings.json``:

   .. code-block:: json

      {
          "cmake.generator": "Ninja",
          "cmake.configureSettings": {
              "CMAKE_BUILD_TYPE": "Debug"
          },
          "cmake.buildDirectory": "${workspaceFolder}/build",
          "cmake.cmakePath": "C:/msys64/mingw64/bin/cmake.exe"
      }

Building the Project
~~~~~~~~~~~~~~~~~~~

1. Open VS Code in the project root directory
2. Press ``Ctrl+Shift+P`` and run "CMake: Configure"
3. Press ``Ctrl+Shift+P`` and run "CMake: Build"

Running Tests
~~~~~~~~~~~~~

In VS Code terminal:

.. code-block:: bash

   cd build
   ctest --output-on-failure

Troubleshooting
~~~~~~~~~~~~~~~

Common issues and solutions:

- **CMake not found**: Ensure MSYS2's MinGW64 bin directory (C:/msys64/mingw64/bin) is in your system PATH
- **Build errors**: Run ``pacman -Syu`` to ensure all packages are up to date
- **CUDA errors**: Install CUDA Toolkit from NVIDIA's website if GPU support is needed 