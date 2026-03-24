#ifndef THREAD_RUNTIME_HPP
#define THREAD_RUNTIME_HPP

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "fun_modern.hpp"

namespace metacsst::app {

inline std::string split_file_path(const std::string& tmp_dir, int index, bool legacy_sub_format) {
  std::ostringstream oss;
  if (legacy_sub_format) {
    if (index >= 10) {
      oss << tmp_dir << "/split_" << index;
    } else {
      oss << tmp_dir << "/split_0" << index;
    }
  } else {
    if (index >= 10) {
      oss << tmp_dir << "/split_" << index;
    } else {
      oss << tmp_dir << "/split_" << std::setw(2) << std::setfill('0') << index;
    }
  }
  return oss.str();
}

inline std::string worker_output_path(const std::string& tmp_dir, int index) {
  std::ostringstream oss;
  oss << tmp_dir << "/out_tmp_" << index << ".txt";
  return oss.str();
}

template <typename ArgType, typename Runner, typename Binder>
inline bool run_split_workers(const std::string& search,
                              int thread_count,
                              const std::string& tmp_dir,
                              bool legacy_sub_format,
                              Binder bind_arg,
                              Runner runner,
                              std::vector<std::string>& output_files) {
  const int number = metacsst::split(search, thread_count, tmp_dir);
  if (number <= 0) {
    return false;
  }

  std::vector<std::thread> threads;
  std::vector<std::unique_ptr<ArgType>> args;
  std::vector<char> worker_success;
  output_files.assign(static_cast<std::size_t>(number), std::string());
  args.reserve(static_cast<std::size_t>(number));
  threads.reserve(static_cast<std::size_t>(number));
  worker_success.assign(static_cast<std::size_t>(number), static_cast<char>(0));

  for (int i = 0; i < number; ++i) {
    const std::string sub_search = split_file_path(tmp_dir, i, legacy_sub_format);
    const std::string out_tmp = worker_output_path(tmp_dir, i);
    const std::size_t index = static_cast<std::size_t>(i);
    output_files[static_cast<std::size_t>(i)] = out_tmp;

    args.push_back(std::make_unique<ArgType>());
    bind_arg(*args.back(), sub_search, out_tmp);
    threads.emplace_back([&, index]() {
      worker_success[index] = runner(*args[index]) ? static_cast<char>(1) : static_cast<char>(0);
    });
  }

  for (auto& t : threads) {
    t.join();
  }

  for (const char ok : worker_success) {
    if (ok == static_cast<char>(0)) {
      return false;
    }
  }
  return true;
}

inline bool merge_files_to_output(const std::vector<std::string>& files, const std::string& output_path, bool append_mode) {
  std::ofstream out;
  if (append_mode) {
    out.open(output_path, std::ios::app);
  } else {
    out.open(output_path);
  }
  if (!out) {
    return false;
  }

  for (const auto& file : files) {
    std::ifstream in(file);
    if (in) {
      out << in.rdbuf();
    }
  }
  return true;
}

inline bool merge_files_to_stream(const std::vector<std::string>& files, std::ostream& out) {
  for (const auto& file : files) {
    std::ifstream in(file);
    if (in) {
      out << in.rdbuf();
    }
  }
  return true;
}

inline void cleanup_tmp_directory(const std::string& tmp_dir) {
  std::filesystem::remove_all(tmp_dir);
}

}

#endif
