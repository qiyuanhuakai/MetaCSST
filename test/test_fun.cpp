#include "test_framework.hpp"
#include "../src/fun_modern.hpp"

#include <vector>
#include <cmath>

using namespace metacsst;

// ==================== String tests ====================

TEST(substr_positive_start) {
    ASSERT_EQ(substr("hello world", 0, 5), "hello");
    ASSERT_EQ(substr("hello world", 6, 5), "world");
    ASSERT_EQ(substr("ATCG", 0, 2), "AT");
    ASSERT_EQ(substr("ATCG", 2, 2), "CG");
    return true;
}

TEST(substr_negative_start) {
    ASSERT_EQ(substr("hello", -1, 1), "o");
    ASSERT_EQ(substr("hello", -2, 2), "ll");
    ASSERT_EQ(substr("ATCG", -1, 1), "G");
    ASSERT_EQ(substr("ATCG", -4, 1), "A");
    return true;
}

TEST(substr_out_of_bounds) {
    ASSERT_EQ(substr("hi", 10, 5), "");
    ASSERT_EQ(substr("hi", -10, 5), "");
    return true;
}

TEST(chomp_basic) {
    std::string s1 = "hello\n";
    ASSERT_EQ(chomp(s1), "hello");

    std::string s2 = "hello";
    ASSERT_EQ(chomp(s2), "hello");

    // Test const version
    ASSERT_EQ(chomp("hello\n"), "hello");
    ASSERT_EQ(chomp("hello"), "hello");
    return true;
}

TEST(chomp_empty) {
    std::string s = "";
    ASSERT_EQ(chomp(s), "");
    ASSERT_EQ(chomp(""), "");
    return true;
}

TEST(judge_contains_gt) {
    ASSERT_EQ(judge(">header"), -1);
    ASSERT_EQ(judge("no gt here"), 0);
    ASSERT_EQ(judge("some>thing"), -1);
    ASSERT_EQ(judge(""), 0);
    return true;
}

// ==================== Sorting tests ====================

TEST(swap_vector) {
    std::vector<float> vec = {1.0f, 2.0f, 3.0f};
    swap(vec, 0, 2);
    ASSERT_EQ(vec[0], 3.0f);
    ASSERT_EQ(vec[2], 1.0f);
    ASSERT_EQ(vec[1], 2.0f);
    return true;
}

TEST(q_sort_vector) {
    std::vector<float> vec = {5.0f, 2.0f, 8.0f, 1.0f, 9.0f};
    q_sort(vec, 0, 4);
    ASSERT_EQ(vec[0], 1.0f);
    ASSERT_EQ(vec[1], 2.0f);
    ASSERT_EQ(vec[2], 5.0f);
    ASSERT_EQ(vec[3], 8.0f);
    ASSERT_EQ(vec[4], 9.0f);
    return true;
}

TEST(q_sort_partial) {
    std::vector<float> vec = {5.0f, 2.0f, 8.0f, 1.0f, 9.0f};
    q_sort(vec, 1, 3);
    ASSERT_EQ(vec[0], 5.0f);  // Unchanged
    ASSERT_EQ(vec[1], 1.0f);
    ASSERT_EQ(vec[2], 2.0f);
    ASSERT_EQ(vec[3], 8.0f);
    ASSERT_EQ(vec[4], 9.0f);  // Unchanged
    return true;
}

// ==================== Cutoff tests ====================

TEST(cutoff_basic) {
    std::vector<float> scores = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    // With ratio 0.5, we want top 50%
    float result = cutoff(scores, 5, 0.5f);
    // Result should be around 3.0 (the middle value after sorting)
    ASSERT_TRUE(result >= 3.0f && result <= 4.0f);
    return true;
}

TEST(cutoff_high_ratio) {
    std::vector<float> scores = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    // With ratio 1.0, we want all values
    float result = cutoff(scores, 5, 1.0f);
    // Should be minimum value
    ASSERT_EQ(result, 1.0f);
    return true;
}

TEST(cutoff_low_ratio) {
    std::vector<float> scores = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    // With ratio 0.0, we want top values
    float result = cutoff(scores, 5, 0.0f);
    // Should be maximum value (or close to it)
    ASSERT_EQ(result, 5.0f);
    return true;
}

TEST(cuttof_alias) {
    // cuttof is an alias for cutoff
    std::vector<float> scores = {1.0f, 2.0f, 3.0f};
    ASSERT_EQ(cuttof(scores, 3, 0.5f), cutoff(scores, 3, 0.5f));
    return true;
}

// ==================== DNA complement tests ====================

TEST(complementary_cases) {
    struct Case { std::string in; std::string out; };
    const std::vector<Case> cases = {
        // AT pairs
        {"AT", "AT"},
        {"TA", "TA"},
        // CG pairs
        {"CG", "CG"},
        {"GC", "GC"},
        // Mixed sequences
        {"ATCG", "CGAT"},
        {"AAA", "TTT"},
        {"TTT", "AAA"},
        {"GGG", "CCC"},
        {"CCC", "GGG"},
        // Case handling (output always uppercase)
        {"atcg", "CGAT"},
        {"AtCg", "CGAT"},
        // With N
        {"ATNCG", "CGNAT"},
        // Empty
        {"", ""}
    };
    
    for (const auto& c : cases) {
        ASSERT_EQ(complementary(c.in), c.out);
    }
    return true;
}

// ==================== Count tests ====================

TEST(count_basic) {
    ASSERT_EQ(count("hello world", 'l'), 3);
    ASSERT_EQ(count("hello world", 'o'), 2);
    ASSERT_EQ(count("hello world", 'x'), 0);
    ASSERT_EQ(count("", 'a'), 0);
    return true;
}

TEST(count_dna) {
    ASSERT_EQ(count("ATCGATCG", 'A'), 2);
    ASSERT_EQ(count("ATCGATCG", 'T'), 2);
    ASSERT_EQ(count("ATCGATCG", 'C'), 2);
    ASSERT_EQ(count("ATCGATCG", 'G'), 2);
    return true;
}

// ==================== Array split tests ====================

TEST(array_split_cases) {
    std::string line = "col1\tcol2\tcol3\tcol4";
    
    // Positive indices
    ASSERT_EQ(array_split(line, '\t', 0), "col1");
    ASSERT_EQ(array_split(line, '\t', 1), "col2");
    ASSERT_EQ(array_split(line, '\t', 2), "col3");
    ASSERT_EQ(array_split(line, '\t', 3), "col4");
    
    // Negative indices
    ASSERT_EQ(array_split(line, '\t', -1), "col3");
    ASSERT_EQ(array_split(line, '\t', -2), "col2");
    ASSERT_EQ(array_split(line, '\t', -3), "col1");
    ASSERT_EQ(array_split(line, '\t', -4), "col1");
    
    // Out of bounds - returns last/first
    std::string short_line = "col1\tcol2";
    ASSERT_EQ(array_split(short_line, '\t', 10), "col2");
    ASSERT_EQ(array_split(short_line, '\t', -10), "col1");
    
    // Single column
    ASSERT_EQ(array_split("only", '\t', 0), "only");
    ASSERT_EQ(array_split("only", '\t', -1), "only");
    
    return true;
}

// ==================== Min/Max tests ====================

TEST(tri_min_max) {
    // Max tests
    ASSERT_EQ(tri_max(1, 2, 3), 3);
    ASSERT_EQ(tri_max(3, 2, 1), 3);
    ASSERT_EQ(tri_max(2, 3, 1), 3);
    ASSERT_EQ(tri_max(5, 5, 5), 5);
    ASSERT_EQ(tri_max(-1, -2, -3), -1);
    
    // Min tests
    ASSERT_EQ(tri_min(1, 2, 3), 1);
    ASSERT_EQ(tri_min(3, 2, 1), 1);
    ASSERT_EQ(tri_min(2, 1, 3), 1);
    ASSERT_EQ(tri_min(5, 5, 5), 5);
    ASSERT_EQ(tri_min(-1, -2, -3), -3);
    
    return true;
}

// ==================== Arg name tests ====================

TEST(arg_name_mapping) {
    // Known mappings
    ASSERT_EQ(arg_name("cov"), "-cov");
    ASSERT_EQ(arg_name("len"), "-len");
    ASSERT_EQ(arg_name("score"), "-score");
    ASSERT_EQ(arg_name("ratio"), "-ratio");
    ASSERT_EQ(arg_name("gap"), "-gap");
    ASSERT_EQ(arg_name("motif"), "-build");
    
    // Unknown mappings
    ASSERT_EQ(arg_name("unknown"), "");
    ASSERT_EQ(arg_name(""), "");
    
    return true;
}

// ==================== Gzip tests ====================

TEST(has_gzip_extension) {
    ASSERT_TRUE(has_gzip_extension("file.gz"));
    ASSERT_TRUE(has_gzip_extension("/path/to/file.gz"));
    ASSERT_FALSE(has_gzip_extension("file.txt"));
    ASSERT_FALSE(has_gzip_extension("file"));
    ASSERT_FALSE(has_gzip_extension("gz"));
    ASSERT_FALSE(has_gzip_extension(""));
    return true;
}

int main() {
    RUN_TESTS();
}
