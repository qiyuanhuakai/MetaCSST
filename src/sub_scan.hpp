#ifndef SUB_SCAN_HPP
#define SUB_SCAN_HPP

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "fasta_runtime.hpp"
#include "ghmm_modern.hpp"
#include "thread_runtime.hpp"

namespace metacsst::subscan {

struct ScanArgs {
  std::string search;
  std::string output;
  HMM_class model;
};

inline void scan_worker(ScanArgs& args) {
  std::ofstream out(args.output);
  if (!out) {
    std::cerr << "Error: Cannot create output file: " << args.output << std::endl;
    return;
  }

  if (args.search.empty()) {
    return;
  }

  const int status = metacsst::app::stream_fasta_sequences(
    args.search,
    true,
    [&](const std::string& name, const std::string& line) {
      auto result = args.model.scanSeq(line);
      if (result->number > 0) {
        out << name << "\n";
        out << line << "\n";
        for (int i = 0; i < result->number; ++i) {
          const std::string match_seq = metacsst::substr(line, result->start[i],
                                                         result->end[i] - result->start[i] + 1);

          if (result->string[i] == 1) {
            out << "Score:" << std::fixed << std::setprecision(2) << result->score[i]
                << "\t+\tmatchSeq(" << result->start[i] << "-" << result->end[i]
                << "):" << match_seq << "\n";
          } else {
            const std::string complementary_seq = metacsst::complementary(match_seq);
            out << "Score:" << std::fixed << std::setprecision(2) << result->score[i]
                << "\t-\tmatchSeq(" << result->start[i] << "-" << result->end[i]
                << "):" << complementary_seq << "\n";
          }
        }
      }
    }
  );

  if (status != 0) {
    return;
  }
}

inline bool run_scan_pipeline(const std::string& search,
                              const std::string& output_path,
                              const std::string& tmp_dir,
                              int thread_count,
                              const HMM_class& model) {
  namespace fs = std::filesystem;

  if (thread_count == 1) {
    const std::string out_tmp = tmp_dir + "/out_tmp.txt";
    ScanArgs arg;
    arg.search = search;
    arg.output = out_tmp;
    arg.model = model;
    scan_worker(arg);

    if (!fs::exists(out_tmp)) {
      std::cerr << "Error: Temporary output file not created: " << out_tmp << std::endl;
      return false;
    }

    fs::copy_file(out_tmp, output_path, fs::copy_options::overwrite_existing);
    return true;
  }

  std::vector<std::string> out_tmp_files;
  const bool worker_ok = metacsst::app::run_split_workers<ScanArgs>(
    search,
    thread_count,
    tmp_dir,
    true,
    [&](ScanArgs& arg, const std::string& split_file, const std::string& worker_output) {
      arg.search = split_file;
      arg.output = worker_output;
      arg.model = model;
    },
    scan_worker,
    out_tmp_files
  );

  if (!worker_ok) {
    std::cerr << "Error: Failed to create split workers." << std::endl;
    return false;
  }

  std::ofstream outfile(output_path);
  if (!outfile) {
    std::cerr << "Error: Cannot create output file: " << output_path << std::endl;
    return false;
  }

  if (!metacsst::app::merge_files_to_stream(out_tmp_files, outfile)) {
    return false;
  }
  return true;
}

}

#endif
