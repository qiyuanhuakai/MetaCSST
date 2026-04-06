#include "test_framework.hpp"
#include "../src/fasta_runtime.hpp"

#include <fstream>
#include <filesystem>
#include <cstdio>

namespace fs = std::filesystem;
using namespace metacsst::app;

// Helper to create temp FASTA files
struct TempFasta {
    std::string path;

    TempFasta(const std::string& filename, const std::string& content) {
        path = (fs::temp_directory_path() / filename).string();
        std::ofstream out(path);
        out << content;
    }

    ~TempFasta() {
        std::remove(path.c_str());
    }
};

// ==================== extract_fasta_name tests ====================

TEST(extract_fasta_name_cases) {
    struct Case {
        std::string header;
        bool stop_at_bracket;
        std::string expected;
    };
    
    const std::vector<Case> cases = {
        // Basic cases
        {">seq1", false, "seq1"},
        {">sequence_name", false, "sequence_name"},
        // With description (stops at space)
        {">seq1 some description", false, "seq1"},
        // Tab is NOT a delimiter (current implementation)
        {">seq1\twith tabs", false, "seq1\twith"},
        // Newline handling
        {">seq1\n", false, "seq1"},
        {">seq1 more\n", false, "seq1"},
        // Bracket handling
        {">seq1[extra]", true, "seq1"},
        {">seq1[extra] desc", true, "seq1"},
        {">seq1[extra]", false, "seq1[extra]"},
        // Edge cases
        {"", false, ""},
        {"not header", false, ""},
        {">", false, ""},
    };
    
    for (const auto& c : cases) {
        ASSERT_EQ(extract_fasta_name(c.header, c.stop_at_bracket), c.expected);
    }
    return true;
}

// ==================== stream_fasta_sequences tests ====================

TEST(stream_single_sequence) {
    TempFasta fa("single.fa",
        ">seq1\n"
        "ATCGATCG\n"
        "ATCG\n"
    );

    int call_count = 0;
    std::string last_name;
    std::string last_seq;

    int result = stream_fasta_sequences(fa.path, false,
        [&](const std::string& name, const std::string& seq) {
            ++call_count;
            last_name = name;
            last_seq = seq;
        }
    );

    ASSERT_EQ(result, 0);
    ASSERT_EQ(call_count, 2);
    ASSERT_EQ(last_name, "seq1");
    ASSERT_EQ(last_seq, "ATCG");
    return true;
}

TEST(stream_multiple_sequences) {
    TempFasta fa("multi.fa",
        ">seq1\n"
        "AAA\n"
        ">seq2\n"
        "TTT\n"
        ">seq3\n"
        "CCC\n"
    );

    std::vector<std::pair<std::string, std::string>> sequences;

    int result = stream_fasta_sequences(fa.path, false,
        [&](const std::string& name, const std::string& seq) {
            sequences.push_back({name, seq});
        }
    );

    ASSERT_EQ(result, 0);
    ASSERT_EQ(sequences.size(), 3);
    ASSERT_EQ(sequences[0].first, "seq1");
    ASSERT_EQ(sequences[0].second, "AAA");
    ASSERT_EQ(sequences[1].first, "seq2");
    ASSERT_EQ(sequences[1].second, "TTT");
    ASSERT_EQ(sequences[2].first, "seq3");
    ASSERT_EQ(sequences[2].second, "CCC");
    return true;
}

TEST(stream_multiline_sequence) {
    TempFasta fa("multiline.fa",
        ">seq1\n"
        "ATCG\n"
        "ATCG\n"
        "ATCG\n"
    );

    int call_count = 0;
    bool all_names_correct = true;

    int result = stream_fasta_sequences(fa.path, false,
        [&](const std::string& name, const std::string& seq) {
            (void)seq;
            ++call_count;
            if (name != "seq1") {
                all_names_correct = false;
            }
        }
    );

    ASSERT_EQ(result, 0);
    ASSERT_EQ(call_count, 3);
    ASSERT_TRUE(all_names_correct);
    return true;
}

TEST(stream_empty_lines) {
    TempFasta fa("emptylines.fa",
        ">seq1\n"
        "\n"
        "ATCG\n"
        "\n"
        ">seq2\n"
        "\n"
        "GGGG\n"
    );

    std::vector<std::pair<std::string, std::string>> sequences;

    int result = stream_fasta_sequences(fa.path, false,
        [&](const std::string& name, const std::string& seq) {
            if (!seq.empty()) {
                sequences.push_back({name, seq});
            }
        }
    );

    ASSERT_EQ(result, 0);
    ASSERT_EQ(sequences.size(), 2);
    ASSERT_EQ(sequences[0].second, "ATCG");
    ASSERT_EQ(sequences[1].second, "GGGG");
    return true;
}

TEST(stream_no_header) {
    TempFasta fa("noheader.fa",
        "ATCG\n"
        "ATCG\n"
    );

    int call_count = 0;

    int result = stream_fasta_sequences(fa.path, false,
        [&](const std::string& name, const std::string& seq) {
            (void)name;
            (void)seq;
            ++call_count;
        }
    );

    ASSERT_EQ(result, 0);
    ASSERT_EQ(call_count, 2);
    return true;
}

TEST(stream_with_bracket_names) {
    TempFasta fa("bracket.fa",
        ">seq1[organism=A] description\n"
        "ATCG\n"
        ">seq2[organism=B]\n"
        "GGGG\n"
    );

    std::vector<std::string> names_with_bracket;
    std::vector<std::string> names_without_bracket;

    stream_fasta_sequences(fa.path, false,
        [&](const std::string& name, const std::string& seq) {
            (void)seq;
            names_with_bracket.push_back(name);
        }
    );

    stream_fasta_sequences(fa.path, true,
        [&](const std::string& name, const std::string& seq) {
            (void)seq;
            names_without_bracket.push_back(name);
        }
    );

    ASSERT_EQ(names_with_bracket.size(), 2);
    ASSERT_EQ(names_without_bracket.size(), 2);

    // With stop_at_bracket=false, name includes [organism=...]
    ASSERT_TRUE(names_with_bracket[0].find("[organism=A]") != std::string::npos);

    // With stop_at_bracket=true, name stops before [
    ASSERT_EQ(names_without_bracket[0], "seq1");
    ASSERT_EQ(names_without_bracket[1], "seq2");
    return true;
}

TEST(stream_nonexistent_file) {
    int result = stream_fasta_sequences("/nonexistent/path/file.fa", false,
        [&](const std::string& name, const std::string& seq) {
            (void)name;
            (void)seq;
        }
    );

    ASSERT_EQ(result, 1);
    return true;
}

TEST(stream_empty_or_header_only) {
    // Empty file
    {
        TempFasta fa("empty.fa", "");
        int call_count = 0;
        int result = stream_fasta_sequences(fa.path, false,
            [&](const std::string& name, const std::string& seq) {
                (void)name;
                (void)seq;
                ++call_count;
            }
        );
        ASSERT_EQ(result, 0);
        ASSERT_EQ(call_count, 0);
    }
    
    // Only header, no sequence
    {
        TempFasta fa("onlyheader.fa", ">seq1\n");
        int call_count = 0;
        int result = stream_fasta_sequences(fa.path, false,
            [&](const std::string& name, const std::string& seq) {
                (void)name;
                (void)seq;
                ++call_count;
            }
        );
        ASSERT_EQ(result, 0);
        ASSERT_EQ(call_count, 0);
    }
    
    return true;
}

int main() {
    RUN_TESTS();
}
