#ifndef CONFIG_MODERN_HPP
#define CONFIG_MODERN_HPP

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace metacsst::config {

using OrderedKeyValues = std::vector<std::pair<std::string, std::string>>;

inline std::string trim(const std::string& s) {
    std::size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])) != 0) {
        start++;
    }
    std::size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])) != 0) {
        end--;
    }
    return s.substr(start, end - start);
}

inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

inline std::string strip_quotes(const std::string& raw) {
    std::string v = trim(raw);
    if (!v.empty() && v.back() == ',') {
        v.pop_back();
        v = trim(v);
    }
    if (v.size() >= 2) {
        if ((v.front() == '"' && v.back() == '"') || (v.front() == '\'' && v.back() == '\'')) {
            return v.substr(1, v.size() - 2);
        }
    }
    return v;
}

inline std::string remove_inline_comment(const std::string& line, char comment_char) {
    bool in_single = false;
    bool in_double = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"' && !in_single) {
            in_double = !in_double;
        } else if (c == '\'' && !in_double) {
            in_single = !in_single;
        } else if (c == comment_char && !in_single && !in_double) {
            return line.substr(0, i);
        }
    }
    return line;
}

inline std::string read_all(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

inline std::string ext_lower(const std::string& path) {
    return to_lower(std::filesystem::path(path).extension().string());
}

inline void ensure_supported_extension(const std::string& path) {
    const std::string ext = ext_lower(path);
    if (ext == ".json" || ext == ".toml" || ext == ".yaml" || ext == ".yml") {
        return;
    }
    throw std::runtime_error("Unsupported config format: " + path + " (only .json/.toml/.yaml/.yml are supported)");
}

inline std::string resolve_path(const std::string& parent_config, const std::string& candidate) {
    std::filesystem::path p(candidate);
    if (p.is_absolute()) {
        return p.lexically_normal().string();
    }
    const auto base = std::filesystem::path(parent_config).parent_path();
    const auto from_config_dir = (base / p).lexically_normal();
    if (std::filesystem::exists(from_config_dir)) {
        return from_config_dir.string();
    }

    const auto config_parent = base.parent_path();
    if (!config_parent.empty()) {
        const auto from_config_parent = (config_parent / p).lexically_normal();
        if (std::filesystem::exists(from_config_parent)) {
            return from_config_parent.string();
        }
    }

    const auto from_cwd = std::filesystem::path(candidate).lexically_normal();
    if (std::filesystem::exists(from_cwd)) {
        return from_cwd.string();
    }

    return from_config_dir.string();
}

inline std::string json_extract_simple_string(const std::string& content, const std::string& key) {
    const std::regex direct("\\\"" + key + "\\\"\\s*:\\s*\\\"([^\\\"]+)\\\"");
    std::smatch m;
    if (std::regex_search(content, m, direct)) {
        return m[1].str();
    }

    const std::regex nested("\\\"" + key + "\\\"\\s*:\\s*\\{[^{}]*\\\"(config|path|file)\\\"\\s*:\\s*\\\"([^\\\"]+)\\\"[^{}]*\\}");
    if (std::regex_search(content, m, nested)) {
        return m[2].str();
    }
    return "";
}

inline std::vector<std::string> json_split_object_array(const std::string& content, std::size_t array_start_pos) {
    std::vector<std::string> objects;
    bool in_string = false;
    bool escaped = false;
    int brace = 0;
    int bracket = 0;
    std::size_t object_start = std::string::npos;

    for (std::size_t i = array_start_pos; i < content.size(); ++i) {
        const char c = content[i];
        if (in_string) {
            if (escaped) {
                escaped = false;
            } else if (c == '\\') {
                escaped = true;
            } else if (c == '"') {
                in_string = false;
            }
            continue;
        }

        if (c == '"') {
            in_string = true;
            continue;
        }
        if (c == '[') {
            bracket++;
            continue;
        }
        if (c == ']') {
            bracket--;
            if (bracket == 0) {
                break;
            }
            continue;
        }
        if (c == '{') {
            if (brace == 0) {
                object_start = i;
            }
            brace++;
            continue;
        }
        if (c == '}') {
            brace--;
            if (brace == 0 && object_start != std::string::npos) {
                objects.push_back(content.substr(object_start, i - object_start + 1));
                object_start = std::string::npos;
            }
        }
    }
    return objects;
}

inline OrderedKeyValues parse_json_object_kv(const std::string& object_content) {
    OrderedKeyValues kvs;
    const std::regex pair_regex("\\\"([A-Za-z_][A-Za-z0-9_]*)\\\"\\s*:\\s*(\\\"([^\\\"]*)\\\"|[-+]?[0-9]*\\.?[0-9]+)");
    auto begin = std::sregex_iterator(object_content.begin(), object_content.end(), pair_regex);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const std::smatch& m = *it;
        std::string value = m[3].matched ? m[3].str() : m[2].str();
        kvs.emplace_back(m[1].str(), value);
    }
    return kvs;
}

inline std::string normalize_dgr_section_name(std::string section) {
    section = trim(section);
    if (section == "TR" || section == "VR" || section == "RT") {
        return section;
    }
    const std::string suffix = ".motifs";
    if (section.size() > suffix.size() && section.compare(section.size() - suffix.size(), suffix.size(), suffix) == 0) {
        section = section.substr(0, section.size() - suffix.size());
    }
    if (section == "TR" || section == "VR" || section == "RT") {
        return section;
    }
    return "";
}

inline std::string json_extract_key_value_block(const std::string& content, const std::string& key) {
    const std::string token = "\"" + key + "\"";
    const std::size_t key_pos = content.find(token);
    if (key_pos == std::string::npos) {
        return "";
    }

    std::size_t pos = content.find(':', key_pos + token.size());
    if (pos == std::string::npos) {
        return "";
    }
    pos++;
    while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos])) != 0) {
        pos++;
    }
    if (pos >= content.size()) {
        return "";
    }

    const char start = content[pos];
    if (start == '{' || start == '[') {
        const char open = start;
        const char close = (open == '{') ? '}' : ']';
        bool in_string = false;
        bool escaped = false;
        int depth = 0;

        for (std::size_t i = pos; i < content.size(); ++i) {
            const char c = content[i];
            if (in_string) {
                if (escaped) {
                    escaped = false;
                } else if (c == '\\') {
                    escaped = true;
                } else if (c == '"') {
                    in_string = false;
                }
                continue;
            }
            if (c == '"') {
                in_string = true;
                continue;
            }
            if (c == open) {
                depth++;
            } else if (c == close) {
                depth--;
                if (depth == 0) {
                    return content.substr(pos, i - pos + 1);
                }
            }
        }
        return "";
    }

    if (start == '"') {
        bool escaped = false;
        for (std::size_t i = pos + 1; i < content.size(); ++i) {
            const char c = content[i];
            if (escaped) {
                escaped = false;
                continue;
            }
            if (c == '\\') {
                escaped = true;
                continue;
            }
            if (c == '"') {
                return content.substr(pos, i - pos + 1);
            }
        }
    }

    std::size_t end = pos;
    while (end < content.size() && content[end] != ',' && content[end] != '}' && content[end] != '\n') {
        end++;
    }
    return trim(content.substr(pos, end - pos));
}

inline std::vector<OrderedKeyValues> parse_json_motif_groups_block(const std::string& block, const std::string& path, const std::string& section) {
    if (block.empty()) {
        throw std::runtime_error("Section " + section + " is empty in " + path);
    }

    std::vector<OrderedKeyValues> groups;
    const std::string t = trim(block);
    if (t.empty()) {
        throw std::runtime_error("Section " + section + " is empty in " + path);
    }

    if (t.front() == '[') {
        for (const auto& object_text : json_split_object_array(t, 0)) {
            OrderedKeyValues kvs = parse_json_object_kv(object_text);
            if (!kvs.empty()) {
                groups.push_back(std::move(kvs));
            }
        }
    } else if (t.front() == '{') {
        std::size_t motifs_pos = t.find("\"motifs\"");
        if (motifs_pos == std::string::npos) {
            throw std::runtime_error("Section " + section + " must contain motifs array in " + path);
        }
        const std::size_t array_pos = t.find('[', motifs_pos);
        if (array_pos == std::string::npos) {
            throw std::runtime_error("Section " + section + " motifs array parse failed in " + path);
        }
        for (const auto& object_text : json_split_object_array(t, array_pos)) {
            OrderedKeyValues kvs = parse_json_object_kv(object_text);
            if (!kvs.empty()) {
                groups.push_back(std::move(kvs));
            }
        }
    } else {
        throw std::runtime_error("Section " + section + " must be object/array in " + path + " (single-file config only)");
    }

    return groups;
}

inline void validate_dgr_motif_groups(const std::unordered_map<std::string, std::vector<OrderedKeyValues>>& result, const std::string& path) {
    for (const char* key : {"TR", "VR", "RT"}) {
        auto it = result.find(key);
        if (it == result.end() || it->second.empty()) {
            throw std::runtime_error("Main config missing motif groups for section: " + std::string(key) + " in " + path);
        }
        for (const auto& group : it->second) {
            bool has_motif = false;
            for (const auto& kv : group) {
                if (kv.first == "motif") {
                    has_motif = true;
                    break;
                }
            }
            if (!has_motif) {
                throw std::runtime_error("Each motif group in section " + std::string(key) + " must contain \"motif\" field in: " + path);
            }
        }
    }
}

inline void resolve_motif_paths(std::vector<OrderedKeyValues>& groups, const std::string& config_path) {
    for (auto& group : groups) {
        for (auto& kv : group) {
            if (kv.first == "motif") {
                kv.second = resolve_path(config_path, kv.second);
            }
        }
    }
}

inline void resolve_dgr_motif_paths(std::unordered_map<std::string, std::vector<OrderedKeyValues>>& result, const std::string& config_path) {
    for (auto& entry : result) {
        resolve_motif_paths(entry.second, config_path);
    }
}

inline std::unordered_map<std::string, std::vector<OrderedKeyValues>> parse_dgr_motif_groups_json(const std::string& path) {
    const std::string content = read_all(path);
    std::unordered_map<std::string, std::vector<OrderedKeyValues>> result;
    for (const char* key : {"TR", "VR", "RT"}) {
        const std::string block = json_extract_key_value_block(content, key);
        if (block.empty()) {
            throw std::runtime_error("Main config missing required section: " + std::string(key) + " in " + path);
        }
        result[key] = parse_json_motif_groups_block(block, path, key);
    }
    validate_dgr_motif_groups(result, path);
    resolve_dgr_motif_paths(result, path);
    return result;
}

inline std::unordered_map<std::string, std::vector<OrderedKeyValues>> parse_dgr_motif_groups_toml(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    std::unordered_map<std::string, std::vector<OrderedKeyValues>> result;
    std::string current_section;
    OrderedKeyValues current_group;

    const auto flush_group = [&]() {
        if (!current_section.empty() && !current_group.empty()) {
            result[current_section].push_back(current_group);
            current_group.clear();
        }
    };

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        line = remove_inline_comment(line, '#');
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }

        if (t.front() == '[' && t.back() == ']') {
            bool is_array_table = t.size() >= 4 && t[1] == '[' && t[t.size() - 2] == ']';
            std::string table_name;
            if (is_array_table) {
                table_name = trim(t.substr(2, t.size() - 4));
            } else {
                table_name = trim(t.substr(1, t.size() - 2));
            }
            const std::string normalized = normalize_dgr_section_name(table_name);
            flush_group();
            current_section = normalized;
            continue;
        }

        const std::size_t eq = t.find('=');
        if (eq == std::string::npos || current_section.empty()) {
            continue;
        }

        const std::string key = trim(t.substr(0, eq));
        const std::string value = strip_quotes(t.substr(eq + 1));
        if (!key.empty() && !value.empty()) {
            current_group.emplace_back(key, value);
        }
    }

    flush_group();
    validate_dgr_motif_groups(result, path);
    resolve_dgr_motif_paths(result, path);
    return result;
}

inline std::unordered_map<std::string, std::vector<OrderedKeyValues>> parse_dgr_motif_groups_yaml(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    std::unordered_map<std::string, std::vector<OrderedKeyValues>> result;
    std::string current_section;
    OrderedKeyValues current_group;

    const auto flush_group = [&]() {
        if (!current_section.empty() && !current_group.empty()) {
            result[current_section].push_back(current_group);
            current_group.clear();
        }
    };

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        line = remove_inline_comment(line, '#');

        int indent = 0;
        while (indent < static_cast<int>(line.size()) && line[static_cast<std::size_t>(indent)] == ' ') {
            indent++;
        }

        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }

        if (indent == 0 && t.back() == ':') {
            const std::string key = trim(t.substr(0, t.size() - 1));
            const std::string normalized = normalize_dgr_section_name(key);
            flush_group();
            current_section = normalized;
            continue;
        }

        if (current_section.empty()) {
            continue;
        }

        if (t.rfind("- ", 0) == 0) {
            flush_group();
            const std::string rest = trim(t.substr(2));
            if (!rest.empty()) {
                const std::size_t colon = rest.find(':');
                if (colon != std::string::npos) {
                    const std::string key = trim(rest.substr(0, colon));
                    const std::string value = strip_quotes(rest.substr(colon + 1));
                    if (!key.empty() && !value.empty()) {
                        current_group.emplace_back(key, value);
                    }
                }
            }
            continue;
        }

        const std::size_t colon = t.find(':');
        if (colon == std::string::npos) {
            continue;
        }
        const std::string key = trim(t.substr(0, colon));
        if (key == "motifs") {
            continue;
        }

        const std::string value = strip_quotes(t.substr(colon + 1));
        if (!key.empty() && !value.empty()) {
            current_group.emplace_back(key, value);
        }
    }

    flush_group();
    validate_dgr_motif_groups(result, path);
    resolve_dgr_motif_paths(result, path);
    return result;
}

inline std::unordered_map<std::string, std::vector<OrderedKeyValues>> parse_dgr_motif_groups(const std::string& path) {
    ensure_supported_extension(path);
    const std::string ext = ext_lower(path);
    if (ext == ".json") {
        return parse_dgr_motif_groups_json(path);
    }
    if (ext == ".toml") {
        return parse_dgr_motif_groups_toml(path);
    }
    return parse_dgr_motif_groups_yaml(path);
}

inline std::vector<OrderedKeyValues> parse_motif_json(const std::string& path) {
    const std::string content = read_all(path);
    std::size_t array_pos = content.find("\"motifs\"");
    if (array_pos != std::string::npos) {
        array_pos = content.find('[', array_pos);
    } else {
        array_pos = content.find('[');
    }

    if (array_pos == std::string::npos) {
        throw std::runtime_error("JSON motif config requires array root or \"motifs\" array: " + path);
    }

    std::vector<OrderedKeyValues> groups;
    for (const auto& object_text : json_split_object_array(content, array_pos)) {
        OrderedKeyValues kvs = parse_json_object_kv(object_text);
        if (!kvs.empty()) {
            groups.push_back(kvs);
        }
    }

    return groups;
}

inline std::vector<OrderedKeyValues> parse_motif_toml(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    std::vector<OrderedKeyValues> groups;
    OrderedKeyValues current;

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        line = remove_inline_comment(line, '#');
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }

        if ((t == "[[motif]]") || (t == "[[motifs]]") || (t == "[motif]") || (t == "[motifs]")) {
            if (!current.empty()) {
                groups.push_back(current);
                current.clear();
            }
            continue;
        }

        const std::size_t eq = t.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        const std::string key = trim(t.substr(0, eq));
        const std::string value = strip_quotes(t.substr(eq + 1));
        if (!key.empty() && !value.empty()) {
            current.emplace_back(key, value);
        }
    }

    if (!current.empty()) {
        groups.push_back(current);
    }
    return groups;
}

inline std::vector<OrderedKeyValues> parse_motif_yaml(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    std::vector<OrderedKeyValues> groups;
    OrderedKeyValues current;
    bool in_motifs = false;

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        line = remove_inline_comment(line, '#');
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }

        if (t == "motifs:") {
            in_motifs = true;
            continue;
        }

        if (t.rfind("- ", 0) == 0) {
            if (!current.empty()) {
                groups.push_back(current);
                current.clear();
            }
            std::string maybe_kv = trim(t.substr(2));
            const std::size_t colon = maybe_kv.find(':');
            if (colon != std::string::npos) {
                std::string key = trim(maybe_kv.substr(0, colon));
                std::string value = strip_quotes(maybe_kv.substr(colon + 1));
                if (!key.empty() && !value.empty()) {
                    current.emplace_back(key, value);
                }
            }
            continue;
        }

        const std::size_t colon = t.find(':');
        if (colon == std::string::npos) {
            continue;
        }

        std::string key = trim(t.substr(0, colon));
        std::string value = strip_quotes(t.substr(colon + 1));

        if (!value.empty() && (in_motifs || !current.empty())) {
            current.emplace_back(key, value);
        }
    }

    if (!current.empty()) {
        groups.push_back(current);
    }
    return groups;
}

inline std::vector<OrderedKeyValues> parse_motif_config(const std::string& path) {
    ensure_supported_extension(path);
    const std::string ext = ext_lower(path);
    std::vector<OrderedKeyValues> groups;

    if (ext == ".json") {
        groups = parse_motif_json(path);
    } else if (ext == ".toml") {
        groups = parse_motif_toml(path);
    } else {
        groups = parse_motif_yaml(path);
    }

    if (groups.empty()) {
        throw std::runtime_error("No motif groups found in config: " + path);
    }

    for (const auto& group : groups) {
        bool has_motif = false;
        for (const auto& kv : group) {
            if (kv.first == "motif") {
                has_motif = true;
                break;
            }
        }
        if (!has_motif) {
            throw std::runtime_error("Each motif group must contain \"motif\" field in: " + path);
        }
    }

    resolve_motif_paths(groups, path);
    return groups;
}

}

#endif
