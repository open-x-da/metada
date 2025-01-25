#include <cuda_runtime.h>
#include <glog/logging.h>
#include <gtest/gtest.h>

class CUDATest : public ::testing::Test {
 protected:
  void SetUp() override { google::InitGoogleLogging("cuda_test"); }

  void TearDown() override { google::ShutdownGoogleLogging(); }

  // Helper function to check CUDA errors
  void checkCudaError(cudaError_t error, const char* msg) {
    ASSERT_EQ(error, cudaSuccess) << msg << ": " << cudaGetErrorString(error);
  }
};

// Test CUDA device properties and availability
TEST_F(CUDATest, DeviceProperties) {
  int deviceCount;
  checkCudaError(cudaGetDeviceCount(&deviceCount),
                 "Failed to get device count");
  ASSERT_GT(deviceCount, 0) << "No CUDA devices found";

  for (int i = 0; i < deviceCount; ++i) {
    cudaDeviceProp prop;
    checkCudaError(cudaGetDeviceProperties(&prop, i),
                   "Failed to get device properties");

    LOG(INFO) << "Device " << i << ": " << prop.name;
    LOG(INFO) << "  Compute capability: " << prop.major << "." << prop.minor;
    LOG(INFO) << "  Total global memory: "
              << prop.totalGlobalMem / (1024 * 1024) << " MB";
    LOG(INFO) << "  Max threads per block: " << prop.maxThreadsPerBlock;

    // Verify minimum requirements
    EXPECT_GE(prop.major, 3) << "Device compute capability too low";
    EXPECT_GT(prop.totalGlobalMem, 1024 * 1024 * 1024ULL)
        << "Insufficient GPU memory";
  }
}

// Test CUDA memory operations
TEST_F(CUDATest, MemoryOperations) {
  const size_t size = 1024 * sizeof(float);
  float* h_data = new float[1024];
  float* d_data = nullptr;

  // Test allocation
  checkCudaError(cudaMalloc(&d_data, size), "Failed to allocate device memory");
  ASSERT_NE(d_data, nullptr);

  // Test memory set
  checkCudaError(cudaMemset(d_data, 0, size), "Failed to set device memory");

  // Test memory copy
  for (int i = 0; i < 1024; ++i) h_data[i] = static_cast<float>(i);
  checkCudaError(cudaMemcpy(d_data, h_data, size, cudaMemcpyHostToDevice),
                 "Failed to copy to device");

  // Verify copy
  float* h_result = new float[1024];
  checkCudaError(cudaMemcpy(h_result, d_data, size, cudaMemcpyDeviceToHost),
                 "Failed to copy from device");

  for (int i = 0; i < 1024; ++i) {
    EXPECT_FLOAT_EQ(h_data[i], h_result[i]);
  }

  // Cleanup
  delete[] h_data;
  delete[] h_result;
  cudaFree(d_data);
}

// Test CUDA stream operations
TEST_F(CUDATest, StreamOperations) {
  cudaStream_t stream;
  checkCudaError(cudaStreamCreate(&stream), "Failed to create CUDA stream");

  // Test stream synchronization
  checkCudaError(cudaStreamSynchronize(stream), "Failed to synchronize stream");

  // Test stream query
  cudaError_t streamStatus = cudaStreamQuery(stream);
  EXPECT_TRUE(streamStatus == cudaSuccess || streamStatus == cudaErrorNotReady);

  checkCudaError(cudaStreamDestroy(stream), "Failed to destroy stream");
}

// Test CUDA event timing
TEST_F(CUDATest, EventTiming) {
  cudaEvent_t start, stop;
  checkCudaError(cudaEventCreate(&start), "Failed to create start event");
  checkCudaError(cudaEventCreate(&stop), "Failed to create stop event");

  // Record events and measure time
  checkCudaError(cudaEventRecord(start), "Failed to record start event");
  cudaDeviceSynchronize();  // Some work would go here
  checkCudaError(cudaEventRecord(stop), "Failed to record stop event");
  checkCudaError(cudaEventSynchronize(stop),
                 "Failed to synchronize stop event");

  float milliseconds = 0;
  checkCudaError(cudaEventElapsedTime(&milliseconds, start, stop),
                 "Failed to calculate elapsed time");
  LOG(INFO) << "Event timing test took " << milliseconds << " ms";

  // Cleanup
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}