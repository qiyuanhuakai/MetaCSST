#include "test_framework.hpp"
#include "../src/config_modern.hpp"

#include <fstream>
#include <filesystem>
#include <cstdio>

namespace fs = std::filesystem;
using namespace metacsst::config;

// Helper function to create temp config files
struct TempConfig {
    std::string path;

    TempConfig(const std::string& filename, const std::string& content) {
        path = (fs::temp_directory_path() / filename).string();
        std::ofstream out(path);
        out << content;
    }

    ~TempConfig() {
        std::remove(path.c_str());
    }
};

TEST(trim_basic) {
    ASSERT_EQ(trim("  hello  "), "hello");
    ASSERT_EQ(trim("hello"), "hello");
    ASSERT_EQ(trim("   hello"), "hello");
    ASSERT_EQ(trim("hello   "), "hello");
    ASSERT_EQ(trim(""), "");
    ASSERT_EQ(trim("   "), "");
    return true;
}

TEST(trim_with_tabs) {
    ASSERT_EQ(trim("\thello\t"), "hello");
    ASSERT_EQ(trim("\t\thello\t\t"), "hello");
    return true;
}

TEST(to_lower_basic) {
    ASSERT_EQ(to_lower("HELLO"), "hello");
    ASSERT_EQ(to_lower("Hello"), "hello");
    ASSERT_EQ(to_lower("hello"), "hello");
    ASSERT_EQ(to_lower("HELLO WORLD"), "hello world");
    ASSERT_EQ(to_lower(""), "");
    return true;
}

TEST(ext_lower_basic) {
    ASSERT_EQ(ext_lower("test.json"), ".json");
    ASSERT_EQ(ext_lower("test.JSON"), ".json");
    ASSERT_EQ(ext_lower("test.JsOn"), ".json");
    ASSERT_EQ(ext_lower("/path/to/file.TOML"), ".toml");
    ASSERT_EQ(ext_lower("config.YAML"), ".yaml");
    ASSERT_EQ(ext_lower("config.YML"), ".yml");
    return true;
}

TEST(ensure_supported_extension_valid) {
    ASSERT_NO_THROW(ensure_supported_extension("test.json"));
    ASSERT_NO_THROW(ensure_supported_extension("test.toml"));
    ASSERT_NO_THROW(ensure_supported_extension("test.yaml"));
    ASSERT_NO_THROW(ensure_supported_extension("test.yml"));
    ASSERT_NO_THROW(ensure_supported_extension("/path/to/config.json"));
    return true;
}

TEST(ensure_supported_extension_invalid) {
    ASSERT_THROW(ensure_supported_extension("test.txt"), std::runtime_error);
    ASSERT_THROW(ensure_supported_extension("test.xml"), std::runtime_error);
    ASSERT_THROW(ensure_supported_extension("test"), std::runtime_error);
    ASSERT_THROW(ensure_supported_extension("test.conf"), std::runtime_error);
    return true;
}

TEST(format_number_string) {
    std::string s1 = format_number_string(3.14);
    ASSERT_TRUE(s1.find("3.14") != std::string::npos);

    std::string s2 = format_number_string(42.0);
    ASSERT_TRUE(s2.find("42") != std::string::npos);

    std::string s3 = format_number_string(0.0);
    ASSERT_TRUE(s3.find("0") != std::string::npos);
    return true;
}

TEST(motif_key_rank_ordering) {
    // Test that motif has highest priority (lowest rank number)
    ASSERT_TRUE(motif_key_rank("motif") < motif_key_rank("cov"));
    ASSERT_TRUE(motif_key_rank("motif") < motif_key_rank("len"));
    ASSERT_TRUE(motif_key_rank("motif") < motif_key_rank("score"));
    ASSERT_TRUE(motif_key_rank("motif") < motif_key_rank("ratio"));
    ASSERT_TRUE(motif_key_rank("motif") < motif_key_rank("gap"));

    // Test order
    ASSERT_TRUE(motif_key_rank("cov") < motif_key_rank("len"));
    ASSERT_TRUE(motif_key_rank("len") < motif_key_rank("score"));
    ASSERT_TRUE(motif_key_rank("score") < motif_key_rank("ratio"));
    ASSERT_TRUE(motif_key_rank("ratio") < motif_key_rank("gap"));

    // Unknown keys get highest rank
    ASSERT_TRUE(motif_key_rank("unknown") >= motif_key_rank("gap"));
    return true;
}

TEST(normalize_group_key_order) {
    OrderedKeyValues group = {
        {"gap", "10"},
        {"motif", "ATCG"},
        {"cov", "0.9"},
        {"score", "5.0"}
    };

    normalize_group_key_order(group);

    // After sorting: motif, cov, score, gap
    ASSERT_EQ(group[0].first, "motif");
    ASSERT_EQ(group[1].first, "cov");
    ASSERT_EQ(group[2].first, "score");
    ASSERT_EQ(group[3].first, "gap");
    return true;
}

TEST(resolve_path_absolute) {
    // Absolute paths should remain unchanged
    std::string abs_path = "/absolute/path/to/file";
    ASSERT_EQ(resolve_path("/config/path.json", abs_path), abs_path);
    return true;
}

TEST(resolve_path_empty) {
    // Empty path should remain empty
    ASSERT_EQ(resolve_path("/config/path.json", ""), "");
    ASSERT_EQ(resolve_path("/config/path.json", "   "), "");
    return true;
}

TEST(parse_json_basic) {
    // Create a simple JSON config with just motifs (without TR/VR/RT sections)
    // Note: parse_motif_json is for motif-only configs
    TempConfig cfg("test_motif.json", R"({
        "motifs": [
            {"motif": "ATCG", "cov": "0.9"}
        ]
    })");

    auto groups = parse_motif_json(cfg.path);
    ASSERT_EQ(groups.size(), 1);
    ASSERT_TRUE(!groups[0].empty());

    // Find motif value
    bool found_motif = false;
    for (const auto& kv : groups[0]) {
        if (kv.first == "motif") {
            ASSERT_EQ(kv.second, "ATCG");
            found_motif = true;
            break;
        }
    }
    ASSERT_TRUE(found_motif);
    return true;
}

TEST(parse_json_array_of_objects) {
    TempConfig cfg("test_array.json", R"([
        {"motif": "AAAA", "cov": "0.8"},
        {"motif": "TTTT", "cov": "0.9"}
    ])");

    auto groups = parse_motif_json(cfg.path);
    ASSERT_EQ(groups.size(), 2);

    bool found_aaaa = false;
    bool found_tttt = false;
    for (const auto& group : groups) {
        for (const auto& kv : group) {
            if (kv.first == "motif" && kv.second == "AAAA") found_aaaa = true;
            if (kv.first == "motif" && kv.second == "TTTT") found_tttt = true;
        }
    }
    ASSERT_TRUE(found_aaaa);
    ASSERT_TRUE(found_tttt);
    return true;
}

TEST(validate_group_has_motif_valid) {
    OrderedKeyValues group = {{"motif", "ATCG"}, {"cov", "0.9"}};
    ASSERT_NO_THROW(validate_group_has_motif(group, "test", "test.json"));
    return true;
}

TEST(validate_group_has_motif_invalid) {
    OrderedKeyValues group = {{"cov", "0.9"}, {"len", "10"}};  // No motif
    ASSERT_THROW(validate_group_has_motif(group, "test", "test.json"), std::runtime_error);
    return true;
}

TEST(validate_group_has_empty_motif) {
    OrderedKeyValues group = {{"motif", ""}, {"cov", "0.9"}};
    ASSERT_THROW(validate_group_has_motif(group, "test", "test.json"), std::runtime_error);
    return true;
}

int main() {
    RUN_TESTS();
}
