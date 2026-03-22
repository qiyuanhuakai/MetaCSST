#ifndef APP_COMMON_HPP
#define APP_COMMON_HPP

#include <filesystem>
#include <exception>
#include <iostream>
#include <string>
#include <string_view>

namespace metacsst::app {

struct CommonOptions {
  int threads = 1;
  std::string config_path;
  std::string input_path;
  std::string output_dir;
  bool show_help = false;
  bool parse_error = false;
};

inline CommonOptions parse_common_options(int argc, char* argv[], const std::string& default_output_dir) {
  CommonOptions options;
  options.output_dir = default_output_dir;

  for (int i = 1; i < argc; ++i) {
    const std::string_view opt(argv[i]);
    if (opt == "-thread" && i + 1 < argc) {
      try {
        options.threads = std::stoi(argv[++i]);
      } catch (const std::exception&) {
        options.parse_error = true;
      }
    } else if (opt == "-build" && i + 1 < argc) {
      options.config_path = argv[++i];
    } else if (opt == "-in" && i + 1 < argc) {
      options.input_path = argv[++i];
    } else if (opt == "-out" && i + 1 < argc) {
      options.output_dir = argv[++i];
    } else if (opt == "-h") {
      options.show_help = true;
    }
  }

  return options;
}

inline void reset_output_directory(const std::string& dir) {
  namespace fs = std::filesystem;
  if (fs::exists(dir)) {
    std::cout << "directory " << dir << " exists,it will be covered!" << std::endl;
    fs::remove_all(dir);
  }
  fs::create_directories(dir);
}

}

#endif
