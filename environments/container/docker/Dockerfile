# Start from NVIDIA CUDA base image
FROM nvidia/cuda:12.1.0-devel-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PATH=/usr/local/cuda/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    ninja-build \
    git \
    python3 \
    python3-pip \
    python3-dev \
    python3-numpy \
    libboost-all-dev \
    libeigen3-dev \
    gfortran \
    doxygen \
    graphviz \
    wget \
    libgtest-dev \
    libgoogle-glog-dev \
    clang-format \
    lcov \
    libyaml-cpp-dev \
    nlohmann-json3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    numpy \
    scipy \
    matplotlib \
    sphinx \
    sphinx_rtd_theme

# Create workspace directory
WORKDIR /opt/metada

# Copy CMake configuration
COPY CMakeLists.txt .
COPY cmake ./cmake

# Copy source code
COPY src ./src
COPY tests ./tests
COPY docs ./docs

# Build the project
RUN cmake -S . -B build \
    && cmake --build build \
    && ctest --test-dir build --output-on-failure

# Set default command
CMD ["/bin/bash"]