#ifndef CONFIG_MODERN_HPP
#define CONFIG_MODERN_HPP

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>
#include <toml++/toml.hpp>
#include <yaml-cpp/yaml.h>

namespace metacsst::config {

using OrderedKeyValues = std::vector<std::pair<std::string, std::string>>;

inline std::string trim(const std::string& s) {
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])) != 0) {
        ++b;
    }
    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) != 0) {
        --e;
    }
    return s.substr(b, e - b);
}

inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

inline std::string ext_lower(const std::string& path) {
    return to_lower(std::filesystem::path(path).extension().string());
}

inline void ensure_supported_extension(const std::string& path) {
    const std::string ext = ext_lower(path);
    if (ext != ".json" && ext != ".toml" && ext != ".yaml" && ext != ".yml") {
        throw std::runtime_error("Unsupported config extension: " + ext + " (expect .json/.toml/.yaml/.yml)");
    }
}

inline std::string read_all(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

inline std::string resolve_path(const std::string& config_path, const std::string& value) {
    const std::string t = trim(value);
    if (t.empty()) {
        return t;
    }
    std::filesystem::path p(t);
    if (p.is_absolute()) {
        return p.lexically_normal().string();
    }
    const std::filesystem::path config_dir = std::filesystem::path(config_path).parent_path();
    const std::filesystem::path from_config_dir = (config_dir / p).lexically_normal();
    if (std::filesystem::exists(from_config_dir)) {
        return from_config_dir.string();
    }

    const std::filesystem::path config_parent = config_dir.parent_path();
    if (!config_parent.empty()) {
        const std::filesystem::path from_config_parent = (config_parent / p).lexically_normal();
        if (std::filesystem::exists(from_config_parent)) {
            return from_config_parent.string();
        }
    }

    const std::filesystem::path from_cwd = std::filesystem::path(t).lexically_normal();
    if (std::filesystem::exists(from_cwd)) {
        return from_cwd.string();
    }

    return from_config_dir.string();
}

inline int motif_key_rank(const std::string& key) {
    static constexpr std::array<const char*, 6> preferred = {"motif", "cov", "len", "score", "ratio", "gap"};
    for (std::size_t i = 0; i < preferred.size(); ++i) {
        if (key == preferred[i]) {
            return static_cast<int>(i);
        }
    }
    return static_cast<int>(preferred.size());
}

inline void normalize_group_key_order(OrderedKeyValues& group) {
    std::stable_sort(group.begin(), group.end(), [](const auto& lhs, const auto& rhs) {
        return motif_key_rank(lhs.first) < motif_key_rank(rhs.first);
    });
}

inline void normalize_motif_groups_order(std::vector<OrderedKeyValues>& groups) {
    for (auto& group : groups) {
        normalize_group_key_order(group);
    }
}

inline void normalize_dgr_groups_order(std::unordered_map<std::string, std::vector<OrderedKeyValues>>& result) {
    for (auto& [key, groups] : result) {
        (void)key;
        normalize_motif_groups_order(groups);
    }
}

inline std::string format_number_string(double value) {
    std::ostringstream oss;
    oss << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return oss.str();
}

inline void validate_group_has_motif(const OrderedKeyValues& group, const std::string& section, const std::string& path) {
    for (const auto& kv : group) {
        if (kv.first == "motif" && !trim(kv.second).empty()) {
            return;
        }
    }
    throw std::runtime_error("Each motif group in section " + section + " must contain non-empty \"motif\" in: " + path);
}

inline void validate_dgr_motif_groups(const std::unordered_map<std::string, std::vector<OrderedKeyValues>>& result, const std::string& path) {
    for (const char* key : {"TR", "VR", "RT"}) {
        const auto it = result.find(key);
        if (it == result.end() || it->second.empty()) {
            throw std::runtime_error("Main config missing motif groups for section: " + std::string(key) + " in " + path);
        }
        for (const auto& group : it->second) {
            validate_group_has_motif(group, key, path);
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
    for (auto& [key, groups] : result) {
        (void)key;
        resolve_motif_paths(groups, config_path);
    }
}

inline std::string json_scalar_to_string(const nlohmann::ordered_json& value, const std::string& path, const std::string& section, const std::string& field) {
    if (value.is_string()) {
        return value.get<std::string>();
    }
    if (value.is_boolean()) {
        return value.get<bool>() ? "true" : "false";
    }
    if (value.is_number_integer()) {
        return std::to_string(value.get<std::int64_t>());
    }
    if (value.is_number_unsigned()) {
        return std::to_string(value.get<std::uint64_t>());
    }
    if (value.is_number_float()) {
        return format_number_string(value.get<double>());
    }
    throw std::runtime_error("Invalid JSON scalar for section " + section + ", field \"" + field + "\" in " + path);
}

inline bool json_object_all_scalar(const nlohmann::ordered_json& obj) {
    if (!obj.is_object()) {
        return false;
    }
    for (auto it = obj.begin(); it != obj.end(); ++it) {
        const auto& v = it.value();
        if (!(v.is_string() || v.is_boolean() || v.is_number())) {
            return false;
        }
    }
    return !obj.empty();
}

inline OrderedKeyValues json_object_to_kvs(const nlohmann::ordered_json& object, const std::string& path, const std::string& section) {
    if (!object.is_object()) {
        throw std::runtime_error("JSON motif group in section " + section + " must be object in " + path);
    }
    OrderedKeyValues kvs;
    kvs.reserve(object.size());
    for (auto it = object.begin(); it != object.end(); ++it) {
        kvs.emplace_back(it.key(), json_scalar_to_string(it.value(), path, section, it.key()));
    }
    return kvs;
}

inline std::vector<OrderedKeyValues> json_node_to_groups(const nlohmann::ordered_json& node, const std::string& path, const std::string& section, bool require_motifs_array) {
    std::vector<OrderedKeyValues> groups;
    if (node.is_array()) {
        groups.reserve(node.size());
        for (const auto& item : node) {
            groups.push_back(json_object_to_kvs(item, path, section));
        }
        return groups;
    }

    if (!node.is_object()) {
        throw std::runtime_error("Section " + section + " must be JSON object/array in " + path);
    }

    if (node.contains("motifs")) {
        const auto& motifs = node.at("motifs");
        if (!motifs.is_array()) {
            throw std::runtime_error("Section " + section + " motifs must be array in " + path);
        }
        groups.reserve(motifs.size());
        for (const auto& item : motifs) {
            groups.push_back(json_object_to_kvs(item, path, section));
        }
        return groups;
    }

    if (!require_motifs_array || json_object_all_scalar(node)) {
        groups.push_back(json_object_to_kvs(node, path, section));
        return groups;
    }

    throw std::runtime_error("Section " + section + " must contain motifs array in " + path);
}

inline std::string toml_scalar_to_string(const toml::node& node, const std::string& path, const std::string& section, const std::string& field) {
    if (const auto* s = node.as_string()) {
        return std::string(s->get());
    }
    if (const auto* i = node.as_integer()) {
        return std::to_string(i->get());
    }
    if (const auto* f = node.as_floating_point()) {
        return format_number_string(f->get());
    }
    if (const auto* b = node.as_boolean()) {
        return b->get() ? "true" : "false";
    }
    throw std::runtime_error("Invalid TOML scalar for section " + section + ", field \"" + field + "\" in " + path);
}

inline OrderedKeyValues toml_table_to_kvs(const toml::table& tbl, const std::string& path, const std::string& section) {
    OrderedKeyValues kvs;
    kvs.reserve(tbl.size());
    for (const auto& [k, v] : tbl) {
        kvs.emplace_back(std::string(k), toml_scalar_to_string(v, path, section, std::string(k)));
    }
    return kvs;
}

inline bool toml_table_all_scalar(const toml::table& tbl) {
    if (tbl.empty()) {
        return false;
    }
    for (const auto& [key, v] : tbl) {
        (void)key;
        if (!(v.is_string() || v.is_integer() || v.is_floating_point() || v.is_boolean())) {
            return false;
        }
    }
    return true;
}

inline std::vector<OrderedKeyValues> toml_node_to_groups(const toml::node& node, const std::string& path, const std::string& section, bool require_motifs_array) {
    std::vector<OrderedKeyValues> groups;

    if (const auto* arr = node.as_array()) {
        groups.reserve(arr->size());
        for (const auto& item : *arr) {
            const auto* tbl = item.as_table();
            if (tbl == nullptr) {
                throw std::runtime_error("TOML motif group in section " + section + " must be table in " + path);
            }
            groups.push_back(toml_table_to_kvs(*tbl, path, section));
        }
        return groups;
    }

    const auto* tbl = node.as_table();
    if (tbl == nullptr) {
        throw std::runtime_error("Section " + section + " must be TOML table/array in " + path);
    }

    if (const auto motifs = (*tbl)["motifs"]) {
        return toml_node_to_groups(*motifs.node(), path, section, false);
    }

    if (!require_motifs_array || toml_table_all_scalar(*tbl)) {
        groups.push_back(toml_table_to_kvs(*tbl, path, section));
        return groups;
    }

    throw std::runtime_error("Section " + section + " must contain motifs array/table in " + path);
}

inline std::string yaml_scalar_to_string(const YAML::Node& node, const std::string& path, const std::string& section, const std::string& field) {
    if (!node.IsScalar()) {
        throw std::runtime_error("Invalid YAML scalar for section " + section + ", field \"" + field + "\" in " + path);
    }
    return node.as<std::string>();
}

inline OrderedKeyValues yaml_map_to_kvs(const YAML::Node& map, const std::string& path, const std::string& section) {
    if (!map.IsMap()) {
        throw std::runtime_error("YAML motif group in section " + section + " must be map in " + path);
    }
    OrderedKeyValues kvs;
    kvs.reserve(map.size());
    for (auto it = map.begin(); it != map.end(); ++it) {
        const std::string key = it->first.as<std::string>();
        kvs.emplace_back(key, yaml_scalar_to_string(it->second, path, section, key));
    }
    return kvs;
}

inline bool yaml_map_all_scalar(const YAML::Node& map) {
    if (!map.IsMap() || map.size() == 0) {
        return false;
    }
    for (auto it = map.begin(); it != map.end(); ++it) {
        if (!it->second.IsScalar()) {
            return false;
        }
    }
    return true;
}

inline std::vector<OrderedKeyValues> yaml_node_to_groups(const YAML::Node& node, const std::string& path, const std::string& section, bool require_motifs_array) {
    std::vector<OrderedKeyValues> groups;
    if (node.IsSequence()) {
        groups.reserve(node.size());
        for (const auto& item : node) {
            groups.push_back(yaml_map_to_kvs(item, path, section));
        }
        return groups;
    }

    if (!node.IsMap()) {
        throw std::runtime_error("Section " + section + " must be YAML map/sequence in " + path);
    }

    const YAML::Node motifs = node["motifs"];
    if (motifs) {
        return yaml_node_to_groups(motifs, path, section, false);
    }

    if (!require_motifs_array || yaml_map_all_scalar(node)) {
        groups.push_back(yaml_map_to_kvs(node, path, section));
        return groups;
    }

    throw std::runtime_error("Section " + section + " must contain motifs sequence/map in " + path);
}

inline std::unordered_map<std::string, std::vector<OrderedKeyValues>> parse_dgr_motif_groups_json(const std::string& path) {
    const auto root = nlohmann::ordered_json::parse(read_all(path));
    if (!root.is_object()) {
        throw std::runtime_error("Main config must be JSON object in " + path);
    }

    std::unordered_map<std::string, std::vector<OrderedKeyValues>> result;
    for (const char* key : {"TR", "VR", "RT"}) {
        const auto it = root.find(key);
        if (it == root.end()) {
            throw std::runtime_error("Main config missing required section: " + std::string(key) + " in " + path);
        }
        result[key] = json_node_to_groups(*it, path, key, true);
    }

    validate_dgr_motif_groups(result, path);
    resolve_dgr_motif_paths(result, path);
    normalize_dgr_groups_order(result);
    return result;
}

inline std::unordered_map<std::string, std::vector<OrderedKeyValues>> parse_dgr_motif_groups_toml(const std::string& path) {
    const auto root = toml::parse_file(path);
    std::unordered_map<std::string, std::vector<OrderedKeyValues>> result;
    for (const char* key : {"TR", "VR", "RT"}) {
        const auto section = root[key];
        if (!section) {
            throw std::runtime_error("Main config missing required section: " + std::string(key) + " in " + path);
        }
        result[key] = toml_node_to_groups(*section.node(), path, key, true);
    }

    validate_dgr_motif_groups(result, path);
    resolve_dgr_motif_paths(result, path);
    normalize_dgr_groups_order(result);
    return result;
}

inline std::unordered_map<std::string, std::vector<OrderedKeyValues>> parse_dgr_motif_groups_yaml(const std::string& path) {
    const auto root = YAML::LoadFile(path);
    if (!root || !root.IsMap()) {
        throw std::runtime_error("Main config must be YAML map in " + path);
    }

    std::unordered_map<std::string, std::vector<OrderedKeyValues>> result;
    for (const char* key : {"TR", "VR", "RT"}) {
        const auto section = root[key];
        if (!section) {
            throw std::runtime_error("Main config missing required section: " + std::string(key) + " in " + path);
        }
        result[key] = yaml_node_to_groups(section, path, key, true);
    }

    validate_dgr_motif_groups(result, path);
    resolve_dgr_motif_paths(result, path);
    normalize_dgr_groups_order(result);
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
    const auto root = nlohmann::ordered_json::parse(read_all(path));
    return json_node_to_groups(root, path, "motifs", false);
}

inline std::vector<OrderedKeyValues> parse_motif_toml(const std::string& path) {
    const auto root = toml::parse_file(path);
    if (const auto motifs = root["motifs"]) {
        return toml_node_to_groups(*motifs.node(), path, "motifs", false);
    }
    if (const auto motif = root["motif"]) {
        return toml_node_to_groups(*motif.node(), path, "motif", false);
    }
    return toml_node_to_groups(root, path, "motif", false);
}

inline std::vector<OrderedKeyValues> parse_motif_yaml(const std::string& path) {
    const auto root = YAML::LoadFile(path);
    if (!root) {
        throw std::runtime_error("YAML motif config is empty: " + path);
    }
    if (const auto motifs = root["motifs"]) {
        return yaml_node_to_groups(motifs, path, "motifs", false);
    }
    return yaml_node_to_groups(root, path, "motif", false);
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
        validate_group_has_motif(group, "motif", path);
    }

    resolve_motif_paths(groups, path);
    normalize_motif_groups_order(groups);
    return groups;
}

}

#endif
