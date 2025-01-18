#include <cuda_runtime.h>
#include "../src/helpers/helper_cuda.h"
#include <gtest/gtest.h>

/**
 * @brief Test fixture for CUDA device query tests
 * 
 * This fixture provides setup and teardown for testing CUDA device properties
 * and capabilities.
 */
class DeviceQueryTest : public ::testing::Test {
protected:
    /** @brief Set up the test environment */
    void SetUp() override {
        // Setup code
    }

    /** @brief Clean up the test environment */
    void TearDown() override {
        // Cleanup code
    }
};

TEST_F(DeviceQueryTest, CheckDeviceCount) {
    int deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
    ASSERT_EQ(error_id, cudaSuccess);
    ASSERT_GT(deviceCount, 0);
}

TEST_F(DeviceQueryTest, CheckDeviceProperties) {
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    
    for (int dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        
        EXPECT_GT(deviceProp.totalGlobalMem, 0);
        EXPECT_GT(deviceProp.multiProcessorCount, 0);
        EXPECT_GT(deviceProp.maxThreadsPerBlock, 0);
    }
} 