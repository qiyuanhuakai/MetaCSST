#ifndef SUB_SCAN_HPP
#define SUB_SCAN_HPP

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "fasta_runtime.hpp"
#include "ghmm_modern.hpp"
#include "scan_pipeline.hpp"
#include "sub_formatter.hpp"

namespace metacsst::subscan {

class sub_scan_strategy : public metacsst::pipeline::ScanStrategy<hmm_class> {
 public:
  using Context = metacsst::pipeline::StrategyContext<hmm_class>;

  bool run_worker(Context& context) const override {
    std::ofstream out(context.output);
    if (!out) {
      std::cerr << "Error: Cannot create output file: " << context.output << std::endl;
      return false;
    }

    if (context.search.empty()) {
      return true;
    }

    const int status = metacsst::app::stream_fasta_sequences(
      context.search,
      true,
      [&](const std::string& name, const std::string& line) {
        auto result = context.model.scan_seq(line);
        if (result->number > 0) {
          metacsst::formatter::write_sub_header(out, name, line);
          for (int i = 0; i < result->number; ++i) {
            const std::string match_seq = metacsst::substr(line, result->start[i],
                                                           result->end[i] - result->start[i] + 1);
            metacsst::formatter::write_sub_score(
              out,
              result->score[i],
              result->string[i],
              result->start[i],
              result->end[i],
              match_seq
            );
          }
        }
      }
    );

    if (status != 0) {
      return false;
    }
    return true;
  }

  bool finalize_single(const std::string& tmp_output,
                       const std::string& final_output) const override {
    namespace fs = std::filesystem;
    fs::copy_file(tmp_output, final_output, fs::copy_options::overwrite_existing);
    return true;
  }

  bool finalize_multi(const std::vector<std::string>& tmp_outputs,
                      const std::string& final_output) const override {
    std::ofstream outfile(final_output);
    if (!outfile) {
      std::cerr << "Error: Cannot create output file: " << final_output << std::endl;
      return false;
    }

    if (!metacsst::app::merge_files_to_stream(tmp_outputs, outfile)) {
      return false;
    }
    return true;
  }

  bool use_legacy_split_naming() const override {
    return true;
  }
};

using SubScanStrategy = sub_scan_strategy;

inline bool run_scan_pipeline(const std::string& search,
                              const std::string& output_path,
                              const std::string& tmp_dir,
                              int thread_count,
                              const hmm_class& model) {
  const sub_scan_strategy strategy;
  return metacsst::pipeline::run_pipeline(search, output_path, tmp_dir, thread_count, model, strategy);
}

}

#endif
