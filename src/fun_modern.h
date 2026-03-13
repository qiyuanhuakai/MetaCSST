#ifndef FUN_MODERN_H
#define FUN_MODERN_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <cctype>
#include <cstring>
#include <sstream>
#include <unordered_map>

namespace metacsst {

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
    using namespace std::filesystem;
    
    int num = 0;
    if (number == 1) return num;
    
    // Count lines in file
    std::ifstream infile(file);
    if (!infile) {
        std::cerr << "Error: Cannot open file " << file << std::endl;
        return 0;
    }
    
    int line = 0;
    std::string tmp;
    while (std::getline(infile, tmp)) {
        line++;
    }
    infile.close();
    
    line /= 2;  // FASTA files have 2 lines per record
    int per = line / number + 1;
    num = line / per;
    if (line % per != 0) {
        num += 1;
    }
    
    // Create output directory if it doesn't exist
    create_directories(dir);
    
    // Read file and split manually (no system() call)
    infile.open(file);
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
    
    while (std::getline(infile, tmp)) {
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
    infile.close();
    
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
    std::cout << "Usage: " << arg << " -build arg.config [Options]\n\n"
              << "Options\n\n"
              << "-build : Config file to build model\n"
              << "[-thread] : Number of threads, [int], default 1\n"
              << "[-in] : Fasta format file, in which patterns are searched, build a HMM only if not given, [string]\n"
              << "[-out] : OUT Directory, [string], default 'sbcsst_out'\n"
              << "[-h] : GHmmMotifScan User Manual\n";
}

// Legacy C-string overloads for backward compatibility during transition
inline char* arg_name(char* name) {
    static thread_local char buffer[20];
    std::string result = arg_name(std::string(name));
    std::strncpy(buffer, result.c_str(), sizeof(buffer) - 1);
    buffer[sizeof(buffer) - 1] = '\0';
    return buffer;
}

inline void usage(char* arg) {
    usage(std::string(arg));
}

inline int judge(char* p) {
    return judge(std::string(p));
}

inline char* chomp(char* s) {
    std::string str(s);
    chomp(str);
    // Note: This modifies the original string in place
    std::size_t len = str.length();
    s[len] = '\0';
    return s;
}

inline char* array_split(char* s, char sep, int col) {
    static thread_local char buffer[constants::D];
    std::string result = array_split(std::string(s), sep, col);
    std::strncpy(buffer, result.c_str(), sizeof(buffer) - 1);
    buffer[sizeof(buffer) - 1] = '\0';
    return buffer;
}

inline char* complementary(char* s) {
    static thread_local char buffer[constants::N];
    std::string result = complementary(std::string(s));
    std::strncpy(buffer, result.c_str(), sizeof(buffer) - 1);
    buffer[sizeof(buffer) - 1] = '\0';
    return buffer;
}

inline int count(char* s, char sep) {
    return count(std::string(s), sep);
}

inline char* substr(char* s, int start, int num) {
    static thread_local char buffer[constants::N];
    std::string result = substr(std::string(s), start, num);
    std::strncpy(buffer, result.c_str(), sizeof(buffer) - 1);
    buffer[sizeof(buffer) - 1] = '\0';
    return buffer;
}

inline int split(char* file, int number, char* dir) {
    return split(std::string(file), number, std::string(dir));
}

// Swapping function for state arrays (used by q_sort_state)
inline void swap_state(int start[], int end[], float score[], int str[], int m, int n) {
    std::swap(start[m], start[n]);
    std::swap(end[m], end[n]);
    std::swap(str[m], str[n]);
    std::swap(score[m], score[n]);
}

// Quick sort for state arrays
inline void q_sort_state(int start[], int end[], float score[], int str[], int left, int right) {
    if (left >= right) return;
    
    // Use middle element as pivot
    swap_state(start, end, score, str, left, (left + right) / 2);
    int last = left;
    
    for (int i = left + 1; i <= right; i++) {
        if (start[i] < start[left]) {
            swap_state(start, end, score, str, ++last, i);
        }
    }
    
    swap_state(start, end, score, str, left, last);
    q_sort_state(start, end, score, str, left, last - 1);
    q_sort_state(start, end, score, str, last + 1, right);
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
using metacsst::swap_state;
using metacsst::q_sort_state;

#endif // FUN_MODERN_H
