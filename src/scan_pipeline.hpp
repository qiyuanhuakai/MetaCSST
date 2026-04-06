#ifndef SCAN_PIPELINE_HPP
#define SCAN_PIPELINE_HPP

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "thread_runtime.hpp"

namespace metacsst::pipeline {

template <typename ModelType>
struct StrategyContext {
  std::string search;
  std::string output;
  ModelType model;
};

template <typename ModelType>
class ScanStrategy {
 public:
  using Context = StrategyContext<ModelType>;
  virtual ~ScanStrategy() = default;
  virtual bool run_worker(Context& context) const = 0;
  virtual bool finalize_single(const std::string& tmp_output,
                               const std::string& final_output) const = 0;
  virtual bool finalize_multi(const std::vector<std::string>& tmp_outputs,
                              const std::string& final_output) const = 0;
  virtual bool use_legacy_split_naming() const = 0;
};

template <typename ModelType>
inline bool run_pipeline(const std::string& search,
                         const std::string& output_path,
                         const std::string& tmp_dir,
                         int thread_count,
                         const ModelType& model,
                         const ScanStrategy<ModelType>& strategy) {
  namespace fs = std::filesystem;
  using Context = StrategyContext<ModelType>;

  if (thread_count == 1) {
    const std::string out_tmp = tmp_dir + "/out_tmp.txt";
    Context context;
    context.search = search;
    context.output = out_tmp;
    context.model = model;

    if (!strategy.run_worker(context)) {
      return false;
    }

    if (!fs::exists(out_tmp)) {
      std::cerr << "Error: Temporary output file not created: " << out_tmp << std::endl;
      return false;
    }
    return strategy.finalize_single(out_tmp, output_path);
  }

  std::vector<std::string> out_tmp_files;
  const bool worker_ok = metacsst::app::run_split_workers<Context>(
    search,
    thread_count,
    tmp_dir,
    strategy.use_legacy_split_naming(),
    [&](Context& context, const std::string& split_file, const std::string& worker_output) {
      context.search = split_file;
      context.output = worker_output;
      context.model = model;
    },
    [&](Context& context) {
      return strategy.run_worker(context);
    },
    out_tmp_files
  );

  if (!worker_ok) {
    std::cerr << "Error: Failed to create split workers." << std::endl;
    return false;
  }
  return strategy.finalize_multi(out_tmp_files, output_path);
}

}

#endif
