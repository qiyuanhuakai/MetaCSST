#ifndef MAIN_SCAN_HPP
#define MAIN_SCAN_HPP

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "config_modern.hpp"
#include "fasta_runtime.hpp"
#include "ghmm_modern.hpp"
#include "thread_runtime.hpp"

namespace metacsst::mainscan {

struct DGRModels {
  HMM_class TR;
  HMM_class VR;
  HMM_class RT;
};

struct ScanArgs {
  std::string search;
  std::string output;
  SCAN scanner;
};

inline DGRModels build_models_from_config(const std::string& config_path) {
  const auto dgr_cfg = metacsst::config::parse_dgr_motif_groups(config_path);
  DGRModels models;
  models.TR.init_groups(dgr_cfg.at("TR"));
  models.VR.init_groups(dgr_cfg.at("VR"));
  models.RT.init_groups(dgr_cfg.at("RT"));
  return models;
}

inline void scan_worker(ScanArgs& args) {
  std::ofstream out(args.output);

  if (!args.search.empty() && out) {
    const int status = metacsst::app::stream_fasta_sequences(
      args.search,
      false,
      [&](const std::string& name, const std::string& seq_str) {
        auto result = args.scanner.scanSeq(seq_str);

        if (result->index == 1) {
          for (int i = 0; i < result->number; ++i) {
            if (result->type[i] != 2) {
              const int start = result->start[i];
              const int end = result->end[i];
              out << name << "\t";
              switch (result->type[i]) {
                case 1: out << "TR\t"; break;
                case 3: out << "RT\t"; break;
              }

              const std::string match_seq = metacsst::substr(seq_str, start, end - start + 1);
              if (result->string[i] == 1) {
                out << std::fixed << std::setprecision(2) << result->score[i] << "\t+\t" << start << "\t" << end << "\t" << match_seq << "\n";
              } else {
                const std::string complementary_seq = metacsst::complementary(match_seq);
                out << std::fixed << std::setprecision(2) << result->score[i] << "\t-\t" << start << "\t" << end << "\t" << complementary_seq << "\n";
              }
            }
          }
          const std::string match_seq = metacsst::substr(seq_str, result->total_start, result->total_end - result->total_start + 1);
          out << name << "\tDGR\t" << std::fixed << std::setprecision(2) << result->total_score << "\t*\t" << result->total_start << "\t" << result->total_end << "\t" << match_seq << "\n";
        }
      }
    );
    if (status != 0) {
      return;
    }
  }
}

inline bool run_scan_pipeline(const std::string& search,
                              const std::string& output_path,
                              const std::string& tmp_dir,
                              int thread_count,
                              const SCAN& scanner) {
  namespace fs = std::filesystem;

  if (thread_count == 1) {
    const std::string out_tmp = tmp_dir + "/out_tmp.txt";
    ScanArgs args;
    args.search = search;
    args.output = out_tmp;
    args.scanner = scanner;

    scan_worker(args);

    if (!fs::exists(out_tmp)) {
      std::cerr << "Error: Temporary output file not created: " << out_tmp << std::endl;
      return false;
    }

    if (!fs::exists(output_path)) {
      fs::copy(out_tmp, output_path, fs::copy_options::overwrite_existing);
      return true;
    }

    const std::vector<std::string> files = {out_tmp};
    if (!metacsst::app::merge_files_to_output(files, output_path, true)) {
      std::cerr << "Error: Failed to append output file: " << output_path << std::endl;
      return false;
    }
    return true;
  }

  std::vector<std::string> out_tmp_files;
  const bool worker_ok = metacsst::app::run_split_workers<ScanArgs>(
    search,
    thread_count,
    tmp_dir,
    false,
    [&](ScanArgs& arg, const std::string& split_file, const std::string& worker_output) {
      arg.search = split_file;
      arg.output = worker_output;
      arg.scanner = scanner;
    },
    scan_worker,
    out_tmp_files
  );

  if (!worker_ok) {
    std::cerr << "Error: Failed to create split workers." << std::endl;
    return false;
  }

  if (!metacsst::app::merge_files_to_output(out_tmp_files, output_path, false)) {
    std::cerr << "Error: Failed to merge worker outputs into: " << output_path << std::endl;
    return false;
  }

  return true;
}

}

#endif
