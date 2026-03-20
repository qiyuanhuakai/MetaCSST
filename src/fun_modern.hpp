#ifndef FUN_MODERN_HPP
#define FUN_MODERN_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <cctype>
#include <cstring>
#include <sstream>
#include <unordered_map>
#include <zlib.h>

namespace metacsst {

inline bool has_gzip_extension(const std::string& path) {
    return path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

inline bool is_gzip_file(const std::string& path) {
    if (has_gzip_extension(path)) {
        return true;
    }

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }

    unsigned char magic[2] = {0, 0};
    in.read(reinterpret_cast<char*>(magic), 2);
    return in.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
}

class LineReader {
public:
    explicit LineReader(const std::string& path)
        : gzip_mode_(is_gzip_file(path)) {
        if (gzip_mode_) {
            gz_file_ = gzopen(path.c_str(), "rb");
        } else {
            file_stream_.open(path);
        }
    }

    ~LineReader() {
        if (gz_file_ != nullptr) {
            gzclose(gz_file_);
        }
    }

    bool is_open() const {
        return gzip_mode_ ? gz_file_ != nullptr : file_stream_.is_open();
    }

    bool getline(std::string& line) {
        line.clear();
        if (!gzip_mode_) {
            return static_cast<bool>(std::getline(file_stream_, line));
        }

        return gzip_getline(line);
    }

private:
    bool gzip_getline(std::string& line) {
        if (gz_file_ == nullptr) {
            return false;
        }

        constexpr int kBufferSize = 8192;
        char buffer[kBufferSize];
        bool has_data = false;

        while (true) {
            char* read_ptr = gzgets(gz_file_, buffer, kBufferSize);
            if (read_ptr == nullptr) {
                return has_data;
            }

            has_data = true;
            const std::size_t len = std::strlen(buffer);
            if (len == 0) {
                if (gzeof(gz_file_)) {
                    return has_data;
                }
                continue;
            }

            if (buffer[len - 1] == '\n') {
                line.append(buffer, len - 1);
                return true;
            }

            line.append(buffer, len);
            if (gzeof(gz_file_)) {
                return true;
            }
        }
    }

    bool gzip_mode_{false};
    std::ifstream file_stream_;
    gzFile gz_file_{nullptr};
};

// Use namespace to avoid global pollution
namespace constants {
    // Max search space number when scanning for sub HMMs
    constexpr std::size_t S = 100000;
    // Max sequence length
    constexpr std::size_t N = 20000000;
    // Max subPattern length
    constexpr std::size_t P = 1000;
    // Max subPattern number
    constexpr std::size_t M = 1000;
    // Max length of directory or file name
    constexpr std::size_t D = 1000;
}

/**
 * Get substring from sequence s
 * @param s Source string
 * @param start Start position (negative = from end)
 * @param num Number of characters to extract
 * @return Substring
 */
inline std::string substr(const std::string& s, int start, int num) {
    num = std::abs(num);
    
    if (start >= 0) {
        if (static_cast<std::size_t>(start) >= s.length()) {
            return "";
        }
        return s.substr(start, static_cast<std::size_t>(num));
    } else {
        // Negative start: count from end
        int from_end = static_cast<int>(s.length()) + start;
        if (from_end < 0) {
            return "";
        }
        std::string result;
        result.reserve(num);
        for (int i = from_end, j = 0; j < num && i >= 0; i--, j++) {
            result.push_back(s[i]);
        }
        return result;
    }
}

/**
 * Remove trailing newline from string (modifies in place)
 * @param s String to modify
 * @return Reference to modified string
 */
inline std::string& chomp(std::string& s) {
    if (!s.empty() && s.back() == '\n') {
        s.pop_back();
    }
    return s;
}

/**
 * Remove trailing newline (returns new string, doesn't modify input)
 * @param s Source string
 * @return String without trailing newline
 */
inline std::string chomp(const std::string& s) {
    std::string result = s;
    return chomp(result);
}

/**
 * Check if string contains '>' character
 * @param p String to check
 * @return -1 if contains '>', 0 otherwise
 */
inline int judge(const std::string& p) {
    return (p.find('>') != std::string::npos) ? -1 : 0;
}

/**
 * Swap two elements in a vector
 * @param vec Vector containing elements
 * @param m First index
 * @param n Second index
 */
inline void swap(std::vector<float>& vec, std::size_t m, std::size_t n) {
    std::swap(vec[m], vec[n]);
}

/**
 * Quick sort using std::sort (wrapper for compatibility)
 * @param vec Vector to sort
 * @param left Left index
 * @param right Right index
 */
inline void q_sort(std::vector<float>& vec, std::size_t left, std::size_t right) {
    if (left >= right) return;
    std::sort(vec.begin() + left, vec.begin() + right + 1);
}

/**
 * Sort using indices (for compatibility with old code that passed pointer to array)
 * @param score Pointer to first element of array
 * @param left Left index
 * @param right Right index
 */
inline void q_sort(float* score, int left, int right) {
    if (left >= right) return;
    std::sort(score + left, score + right + 1);
}

/**
 * Get cutoff value from sorted scores
 * @param score Vector of scores
 * @param number Number of elements
 * @param ratio Ratio between 0 and 1
 * @return Cutoff value
 */
inline float cuttof(const std::vector<float>& score, int number, float ratio) {
    std::vector<float> sorted(score.begin(), score.begin() + number);
    std::sort(sorted.begin(), sorted.end());
    int i = static_cast<int>(std::ceil((1 - ratio) * number));
    i = std::clamp(i, 0, number - 1);
    return sorted[i];
}

/**
 * Get cutoff value from array (legacy interface)
 * @param score Pointer to score array
 * @param number Number of elements
 * @param ratio Ratio between 0 and 1
 * @return Cutoff value
 */
inline float cuttof(float** score, int number, float ratio) {
    q_sort(*score, 0, number);
    int i = static_cast<int>(std::ceil((1 - ratio) * number));
    i = std::clamp(i, 0, number);
    return (*score)[i];
}

/**
 * Split file into parts for multi-threading
 * @param file Path to input file
 * @param number Number of threads
 * @param dir Output directory
 * @return Number of split files created
 */
inline int split(const std::string& file, int number, const std::string& dir) {
    int num = 0;
    if (number == 1) return num;
    
    // Count lines in file
    LineReader infile(file);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open file " << file << std::endl;
        return 0;
    }
    
    int line = 0;
    std::string tmp;
    while (infile.getline(tmp)) {
        line++;
    }
    
    line /= 2;  // FASTA files have 2 lines per record
    int per = line / number + 1;
    num = line / per;
    if (line % per != 0) {
        num += 1;
    }
    
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(dir);
    
    // Read file and split manually (no system() call)
    LineReader split_reader(file);
    if (!split_reader.is_open()) {
        std::cerr << "Error: Cannot open file " << file << std::endl;
        return 0;
    }

    int file_idx = 0;
    int lines_written = 0;
    int lines_per_file = per * 2;
    
    std::ofstream outfile;
    auto open_new_file = [&]() {
        std::ostringstream oss;
        oss << dir << "/split_" << std::setw(2) << std::setfill('0') << file_idx;
        outfile.open(oss.str());
        if (!outfile) {
            std::cerr << "Error: Cannot create output file " << oss.str() << std::endl;
        }
    };
    
    open_new_file();
    
    while (split_reader.getline(tmp)) {
        outfile << tmp << '\n';
        lines_written++;
        
        if (lines_written >= lines_per_file) {
            outfile.close();
            file_idx++;
            lines_written = 0;
            open_new_file();
        }
    }
    
    outfile.close();
    
    return num;
}

/**
 * Get DNA reverse complement
 * @param s DNA sequence
 * @return Reverse complement string
 */
inline std::string complementary(const std::string& s) {
    std::string result;
    result.reserve(s.length());

    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        char c = std::toupper(*it);
        switch (c) {
            case 'A': result.push_back('T'); break;
            case 'T': result.push_back('A'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'N': result.push_back('N'); break;
            default:  result.push_back(c); break;
        }
    }
    return result;
}

/**
 * Count occurrences of character in string
 * @param s String to search
 * @param sep Character to count
 * @return Number of occurrences
 */
inline int count(const std::string& s, char sep) {
    return static_cast<int>(std::count(s.begin(), s.end(), sep));
}

/**
 * Split string by separator and return specified column
 * @param s String to split
 * @param sep Separator character
 * @param col Column index (negative = from end)
 * @return Requested column value
 */
inline std::string array_split(const std::string& s, char sep, int col) {
    std::vector<std::string> columns;
    std::stringstream ss(s);
    std::string item;
    
    while (std::getline(ss, item, sep)) {
        columns.push_back(item);
    }
    
    if (columns.empty()) {
        return "";
    }
    
    std::size_t column_count = columns.size() - 1;  // Original code behavior
    
    if (col >= 0) {
        if (static_cast<std::size_t>(col) > column_count) {
            return columns.back();
        } else {
            return columns[col];
        }
    } else {
        int idx = static_cast<int>(column_count) + col;
        if (idx < 0) {
            return columns.front();
        } else {
            return columns[idx];
        }
    }
}

/**
 * Return maximum of three integers
 */
inline int tri_max(int a, int b, int c) {
    return std::max({a, b, c});
}

/**
 * Return minimum of three integers
 */
inline int tri_min(int a, int b, int c) {
    return std::min({a, b, c});
}

/**
 * Get argument name mapping
 * @param name Parameter name
 * @return Command-line argument string
 */
inline std::string arg_name(const std::string& name) {
    static const std::unordered_map<std::string, std::string> arg_map = {
        {"cov", "-cov"},
        {"len", "-len"},
        {"score", "-score"},
        {"ratio", "-ratio"},
        {"gap", "-gap"},
        {"motif", "-build"}
    };
    
    auto it = arg_map.find(name);
    if (it != arg_map.end()) {
        return it->second;
    }
    return "";
}

/**
 * Print usage information
 * @param arg Program name
 */
inline void usage(const std::string& arg) {
    std::cout << "Usage: " << arg << " -build config.json|config.toml|config.yaml|config.yml [Options]\n\n"
              << "Options\n\n"
              << "-build : Config file to build model\n"
              << "[-thread] : Number of threads, [int], default 1\n"
              << "[-in] : Fasta format file, in which patterns are searched, build a HMM only if not given, [string]\n"
              << "[-out] : OUT Directory, [string], default 'out_metacsst'\n"
              << "[-h] : GHmmMotifScan User Manual\n";
}

} // namespace metacsst

// Pull constants into global namespace for backward compatibility
using metacsst::constants::S;
using metacsst::constants::N;
using metacsst::constants::P;
using metacsst::constants::M;
using metacsst::constants::D;

// Pull functions into global namespace for backward compatibility
using metacsst::substr;
using metacsst::chomp;
using metacsst::judge;
using metacsst::has_gzip_extension;
using metacsst::is_gzip_file;
using metacsst::LineReader;
using metacsst::swap;
using metacsst::q_sort;
using metacsst::cuttof;
using metacsst::split;
using metacsst::complementary;
using metacsst::count;
using metacsst::array_split;
using metacsst::tri_max;
using metacsst::tri_min;
using metacsst::arg_name;
using metacsst::usage;

#endif
