#include "test_framework.hpp"
#include "../src/app_common.hpp"

#include <filesystem>
#include <cstdio>
#include <fstream>

namespace fs = std::filesystem;
using namespace metacsst::app;

// ==================== parse_common_options tests ====================

TEST(parse_no_args) {
    // Test with minimal args (program name only)
    const char* args[] = {"program"};
    auto opts = parse_common_options(1, const_cast<char**>(args), "default_out");

    ASSERT_EQ(opts.threads, 1);  // Default
    ASSERT_EQ(opts.config_path, "");
    ASSERT_EQ(opts.input_path, "");
    ASSERT_EQ(opts.output_dir, "default_out");
    ASSERT_FALSE(opts.show_help);
    ASSERT_FALSE(opts.parse_error);
    return true;
}

TEST(parse_thread_option) {
    const char* args[] = {"program", "-thread", "4"};
    auto opts = parse_common_options(3, const_cast<char**>(args), "out");

    ASSERT_EQ(opts.threads, 4);
    return true;
}

TEST(parse_thread_invalid) {
    const char* args[] = {"program", "-thread", "abc"};
    auto opts = parse_common_options(3, const_cast<char**>(args), "out");

    ASSERT_TRUE(opts.parse_error);
    return true;
}

TEST(parse_thread_negative) {
    const char* args[] = {"program", "-thread", "-1"};
    auto opts = parse_common_options(3, const_cast<char**>(args), "out");

    ASSERT_EQ(opts.threads, -1);  // Parsed as -1, validation should happen later
    return true;
}

TEST(parse_build_option) {
    const char* args[] = {"program", "-build", "config.json"};
    auto opts = parse_common_options(3, const_cast<char**>(args), "out");

    ASSERT_EQ(opts.config_path, "config.json");
    return true;
}

TEST(parse_input_option) {
    const char* args[] = {"program", "-in", "input.fa"};
    auto opts = parse_common_options(3, const_cast<char**>(args), "out");

    ASSERT_EQ(opts.input_path, "input.fa");
    return true;
}

TEST(parse_output_option) {
    const char* args[] = {"program", "-out", "custom_output"};
    auto opts = parse_common_options(3, const_cast<char**>(args), "default");

    ASSERT_EQ(opts.output_dir, "custom_output");
    return true;
}

TEST(parse_help_option) {
    const char* args[] = {"program", "-h"};
    auto opts = parse_common_options(2, const_cast<char**>(args), "out");

    ASSERT_TRUE(opts.show_help);
    return true;
}

TEST(parse_multiple_options) {
    const char* args[] = {
        "program",
        "-build", "config.json",
        "-in", "input.fa",
        "-out", "results",
        "-thread", "8"
    };
    auto opts = parse_common_options(9, const_cast<char**>(args), "default");

    ASSERT_EQ(opts.config_path, "config.json");
    ASSERT_EQ(opts.input_path, "input.fa");
    ASSERT_EQ(opts.output_dir, "results");
    ASSERT_EQ(opts.threads, 8);
    ASSERT_FALSE(opts.show_help);
    ASSERT_FALSE(opts.parse_error);
    return true;
}

TEST(parse_option_without_value) {
    // Test option at end without value
    const char* args[] = {"program", "-build"};
    auto opts = parse_common_options(2, const_cast<char**>(args), "out");

    // Should not crash, should leave config_path empty
    ASSERT_EQ(opts.config_path, "");
    return true;
}

TEST(parse_thread_without_value) {
    const char* args[] = {"program", "-thread"};
    auto opts = parse_common_options(2, const_cast<char**>(args), "out");

    ASSERT_EQ(opts.threads, 1);  // Should keep default
    return true;
}

// ==================== reset_output_directory tests ====================

struct TempDirCleanup {
    std::string path;
    ~TempDirCleanup() {
        if (fs::exists(path)) {
            fs::remove_all(path);
        }
    }
};

TEST(reset_creates_directory) {
    std::string test_dir = (fs::temp_directory_path() / "test_reset_create").string();
    TempDirCleanup cleanup{test_dir};

    // Remove if exists
    if (fs::exists(test_dir)) {
        fs::remove_all(test_dir);
    }

    ASSERT_FALSE(fs::exists(test_dir));

    reset_output_directory(test_dir);

    ASSERT_TRUE(fs::exists(test_dir));
    ASSERT_TRUE(fs::is_directory(test_dir));
    return true;
}

TEST(reset_clears_existing_directory) {
    std::string test_dir = (fs::temp_directory_path() / "test_reset_clear").string();
    TempDirCleanup cleanup{test_dir};

    // Create directory with a file
    fs::create_directories(test_dir);
    std::ofstream((test_dir + "/existing_file.txt")) << "content";

    ASSERT_TRUE(fs::exists(test_dir + "/existing_file.txt"));

    reset_output_directory(test_dir);

    ASSERT_TRUE(fs::exists(test_dir));
    ASSERT_FALSE(fs::exists(test_dir + "/existing_file.txt"));
    return true;
}

TEST(reset_nested_directory) {
    std::string test_dir = (fs::temp_directory_path() / "test_reset" / "nested" / "dir").string();
    TempDirCleanup cleanup{(fs::temp_directory_path() / "test_reset").string()};

    if (fs::exists(fs::temp_directory_path() / "test_reset")) {
        fs::remove_all(fs::temp_directory_path() / "test_reset");
    }

    reset_output_directory(test_dir);

    ASSERT_TRUE(fs::exists(test_dir));
    return true;
}

// ==================== CommonOptions struct tests ====================

TEST(common_options_defaults) {
    common_options opts;

    ASSERT_EQ(opts.threads, 1);
    ASSERT_EQ(opts.config_path, "");
    ASSERT_EQ(opts.input_path, "");
    ASSERT_EQ(opts.output_dir, "");
    ASSERT_FALSE(opts.show_help);
    ASSERT_FALSE(opts.parse_error);
    return true;
}

int main() {
    RUN_TESTS();
}
