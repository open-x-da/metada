Windows Installation Guide
==========================

Setting up Development Environment for METADA
---------------------------------------------

Prerequisites
~~~~~~~~~~~~~

Before starting, ensure you have the following installed:

1. **NVIDIA GPU** with CUDA support (for GPU acceleration)
2. **CUDA Toolkit** (latest version recommended)
   - Download from `NVIDIA CUDA Toolkit <https://developer.nvidia.com/cuda-downloads>`_
   - Choose the Windows installer that matches your system
3. **Visual Studio Build Tools** (required for compilation)
   - Download from `Visual Studio Downloads <https://visualstudio.microsoft.com/downloads/>`_
   - Select "Desktop development with C++" workload during installation

Installing Visual Studio Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download Visual Studio Code from the `official website <https://code.visualstudio.com/>`_
2. Run the installer and follow the installation wizard
3. Launch VS Code after installation

Required Extensions
~~~~~~~~~~~~~~~~~~~

Install the following VS Code extensions:

1. **C/C++** (Microsoft)
   - Provides C++ language support
   - Install from VS Code marketplace or using command:
     
     ``code --install-extension ms-vscode.cpptools``

2. **CMake Tools** (Microsoft)
   - Provides CMake integration
   - Install from VS Code marketplace or using command:
     
     ``code --install-extension ms-vscode.cmake-tools``

Configuring the Build Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Set up CMake configuration**:
   
   Create or modify ``.vscode/settings.json``:

   .. code-block:: json

      {
          "cmake.configureSettings": {
              "CMAKE_CUDA_ARCHITECTURES": "75",
              "CMAKE_CUDA_COMPILER": "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.1/bin/nvcc.exe"
          },
          "cmake.generator": "Visual Studio 17 2022",
          "cmake.buildDirectory": "${workspaceFolder}/build"
      }

2. **Configure IntelliSense**:
   
   Create or modify ``.vscode/c_cpp_properties.json``:

   .. literalinclude:: ../../../.vscode/c_cpp_properties.json
      :language: json
      :caption: .vscode/c_cpp_properties.json

Building METADA
~~~~~~~~~~~~~~~

1. Open VS Code in the METADA root directory
2. Press ``Ctrl+Shift+P`` and type "CMake: Configure"
3. Select your compiler (Visual Studio 2022)
4. Press ``Ctrl+Shift+P`` and type "CMake: Build"

Running Tests
~~~~~~~~~~~~~

1. Build the project as described above
2. Navigate to the build directory
3. Run ``ctest`` or use VS Code's test explorer

Troubleshooting
~~~~~~~~~~~~~~~

Common issues and solutions:

- **CMake configuration fails**: Verify CUDA path in settings.json
- **Build errors**: Ensure Visual Studio Build Tools are properly installed
- **GPU not detected**: Check NVIDIA driver installation and GPU compatibility 