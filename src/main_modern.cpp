#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

#include "ghmm_modern.hpp"
#include "fun_modern.hpp"
#include "app_common.hpp"
#include "main_scan.hpp"

int main(int argc,char* argv[]){
  if(argc < 3){
    metacsst::usage(argv[0]);
    return 0;
  }

  const auto options = metacsst::app::parse_common_options(argc, argv, "out_metacsst");
  if (options.parse_error) {
    metacsst::usage(argv[0]);
    return 1;
  }
  if (options.show_help) {
    metacsst::usage(argv[0]);
    return 0;
  }

  if(options.threads <= 0) {
    metacsst::usage(argv[0]);
    return 1;
  }

  if(options.config_path.empty()) {
    metacsst::usage(argv[0]);
    return 0;
  }

  const std::string& search = options.input_path;
  const std::string& dir = options.output_dir;
  metacsst::app::reset_output_directory(dir);

  const std::string tmp = dir + "/tmp";
  const std::string out = dir + "/raw.gtf";

  metacsst::mainscan::DGRModels models;
  try {
    models = metacsst::mainscan::build_models_from_config(options.config_path);
  } catch (const std::exception& ex) {
    std::cerr << "Config error: " << ex.what() << std::endl;
    return 1;
  }

  SCAN dgrScan;
  dgrScan.init(models.TR, models.VR, models.RT, 10000);
  dgrScan.print(dir);

  if(search.empty()) {
    return 0;
  }

  fs::create_directories(tmp);

  if(!metacsst::mainscan::run_scan_pipeline(search, out, tmp, options.threads, dgrScan)) {
    metacsst::app::cleanup_tmp_directory(tmp);
    return 1;
  }

  metacsst::app::cleanup_tmp_directory(tmp);
  return 0;
}
