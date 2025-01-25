#include <gtest/gtest.h>
#include <glog/logging.h>
#include <fstream>
#include <string>
#include <filesystem>

class GLogTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize glog with test name
        google::InitGoogleLogging("glog_test");
        
        // Set log directory to build directory
        log_dir = std::filesystem::temp_directory_path() / "glog_test";
        std::filesystem::create_directories(log_dir);
        FLAGS_log_dir = log_dir.string();
        
        // Enable logging to stderr for all levels
        FLAGS_logtostderr = false;
        FLAGS_stderrthreshold = 0;
    }

    void TearDown() override {
        // Clean up log files
        google::ShutdownGoogleLogging();
        std::filesystem::remove_all(log_dir);
    }

    std::filesystem::path log_dir;
};

// Helper function to find log file by severity
bool find_log_file(const std::filesystem::path& dir, const std::string& severity) {
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().string().find(severity) != std::string::npos) {
            return true;
        }
    }
    return false;
}

// Helper function to read log file content
std::string read_log_file(const std::filesystem::path& dir, const std::string& severity) {
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().string().find(severity) != std::string::npos) {
            std::ifstream file(entry.path());
            return std::string((std::istreambuf_iterator<char>(file)),
                              std::istreambuf_iterator<char>());
        }
    }
    return "";
}

// Test basic logging levels
TEST_F(GLogTest, LogLevels) {
    LOG(INFO) << "Info message";
    LOG(WARNING) << "Warning message";
    LOG(ERROR) << "Error message";
    
    // Flush all severity levels
    google::FlushLogFiles(google::GLOG_INFO);
    google::FlushLogFiles(google::GLOG_WARNING);
    google::FlushLogFiles(google::GLOG_ERROR);
    google::FlushLogFiles(google::GLOG_FATAL);

    // Verify log files exist with correct pattern
    EXPECT_TRUE(find_log_file(log_dir, "INFO"));
    EXPECT_TRUE(find_log_file(log_dir, "WARNING"));
    EXPECT_TRUE(find_log_file(log_dir, "ERROR"));
}

// Test conditional logging
TEST_F(GLogTest, ConditionalLogging) {
    const bool should_log = true;
    const bool should_not_log = false;
    
    LOG_IF(INFO, should_log) << "This should be logged";
    LOG_IF(INFO, should_not_log) << "This should not be logged";
    
    // Flush the logs to ensure they're written
    google::FlushLogFiles(google::GLOG_INFO);

    std::string content = read_log_file(log_dir, "INFO");
    
    EXPECT_TRUE(content.find("This should be logged") != std::string::npos);
    EXPECT_TRUE(content.find("This should not be logged") == std::string::npos);
}

// Test CHECK macros
TEST_F(GLogTest, CheckMacros) {
    const int x = 5;
    CHECK_EQ(x, 5) << "x should equal 5";
    CHECK_NE(x, 0) << "x should not equal 0";
    CHECK_GT(x, 0) << "x should be positive";
    
    EXPECT_DEATH({
        CHECK_EQ(x, 0) << "This should cause fatal error";
    }, "Check failed");
}

// Test VLOG
TEST_F(GLogTest, VerboseLogging) {
    FLAGS_v = 1;  // Set verbosity level
    
    VLOG(0) << "VLOG 0 message";  // Should appear
    VLOG(1) << "VLOG 1 message";  // Should appear
    VLOG(2) << "VLOG 2 message";  // Should not appear
    
    // Flush the logs to ensure they're written
    google::FlushLogFiles(google::GLOG_INFO);
    
    std::string content = read_log_file(log_dir, "INFO");
    
    EXPECT_TRUE(content.find("VLOG 0 message") != std::string::npos);
    EXPECT_TRUE(content.find("VLOG 1 message") != std::string::npos);
    EXPECT_TRUE(content.find("VLOG 2 message") == std::string::npos);
} 