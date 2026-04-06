#ifndef MAIN_SCAN_HPP
#define MAIN_SCAN_HPP

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "config_modern.hpp"
#include "fasta_runtime.hpp"
#include "ghmm_modern.hpp"
#include "main_formatter.hpp"
#include "scan_pipeline.hpp"

namespace metacsst::mainscan {

struct dgr_models {
  hmm_class TR;
  hmm_class VR;
  hmm_class RT;
};

inline dgr_models build_models_from_config(const std::string& config_path) {
  const auto dgr_cfg = metacsst::config::parse_dgr_motif_groups(config_path);
  dgr_models models;
  models.TR.init_groups(dgr_cfg.at("TR"));
  models.VR.init_groups(dgr_cfg.at("VR"));
  models.RT.init_groups(dgr_cfg.at("RT"));
  return models;
}

class main_scan_strategy : public metacsst::pipeline::ScanStrategy<scan_model> {
 public:
  using Context = metacsst::pipeline::StrategyContext<scan_model>;

  bool run_worker(Context& context) const override {
    std::ofstream out(context.output);
    if (!out) {
      std::cerr << "Error: Cannot create output file: " << context.output << std::endl;
      return false;
    }

    if (!context.search.empty()) {
      const int status = metacsst::app::stream_fasta_sequences(
        context.search,
        false,
        [&](const std::string& name, const std::string& seq_str) {
          auto result = context.model.scan_seq(seq_str);

          if (result->index == 1) {
            for (int i = 0; i < result->number; ++i) {
              if (result->type[i] != 2) {
                const int start = result->start[i];
                const int end = result->end[i];
                const std::string match_seq = metacsst::substr(seq_str, start, end - start + 1);
                metacsst::formatter::write_main_feature(
                  out,
                  name,
                  result->type[i],
                  result->score[i],
                  result->string[i],
                  start,
                  end,
                  match_seq
                );
              }
            }
            const std::string match_seq = metacsst::substr(seq_str, result->total_start, result->total_end - result->total_start + 1);
            metacsst::formatter::write_main_dgr(
              out,
              name,
              result->total_score,
              result->total_start,
              result->total_end,
              match_seq
            );
          }
        }
      );
      if (status != 0) {
        return false;
      }
    }

    return true;
  }

  bool finalize_single(const std::string& tmp_output,
                       const std::string& final_output) const override {
    namespace fs = std::filesystem;
    if (!fs::exists(final_output)) {
      fs::copy(tmp_output, final_output, fs::copy_options::overwrite_existing);
      return true;
    }

    const std::vector<std::string> files = {tmp_output};
    if (!metacsst::app::merge_files_to_output(files, final_output, true)) {
      std::cerr << "Error: Failed to append output file: " << final_output << std::endl;
      return false;
    }
    return true;
  }

  bool finalize_multi(const std::vector<std::string>& tmp_outputs,
                      const std::string& final_output) const override {
    if (!metacsst::app::merge_files_to_output(tmp_outputs, final_output, false)) {
      std::cerr << "Error: Failed to merge worker outputs into: " << final_output << std::endl;
      return false;
    }
    return true;
  }

  bool use_legacy_split_naming() const override {
    return false;
  }
};

using MainScanStrategy = main_scan_strategy;

inline bool run_scan_pipeline(const std::string& search,
                              const std::string& output_path,
                              const std::string& tmp_dir,
                              int thread_count,
                              const scan_model& scanner) {
  const main_scan_strategy strategy;
  return metacsst::pipeline::run_pipeline(search, output_path, tmp_dir, thread_count, scanner, strategy);
}

using DGRModels = dgr_models;

}

#endif
