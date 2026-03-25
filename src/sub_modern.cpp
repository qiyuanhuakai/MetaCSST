#include <filesystem>
#include <iostream>
#include <string>

#include "ghmm_modern.hpp"
#include "app_common.hpp"
#include "sub_scan.hpp"

namespace fs = std::filesystem;

int main(int argc,char* argv[]){
  if(argc < 3){
    metacsst::usage(argv[0]);
    return 0;
  }

  const auto options = metacsst::app::parse_common_options(argc, argv, "sbcsst_out");
  if (options.parse_error) {
    metacsst::usage(argv[0]);
    return 1;
  }
  if (options.show_help) {
    metacsst::usage(argv[0]);
    return 0;
  }

  if(options.threads <= 0){
    metacsst::usage(argv[0]);
    return 1;
  }

  if(options.config_path.empty()){
    metacsst::usage(argv[0]);
    return 0;
  }

  const std::string& search = options.input_path;
  const std::string& dir = options.output_dir;
  metacsst::app::reset_output_directory(dir);

  const std::string tmp = dir + "/tmp";
  const std::string out = dir + "/out.txt";

  hmm_class hmm_class;
  try {
    hmm_class.init(options.config_path);
  } catch (const std::exception& ex) {
    std::cerr << "Config error: " << ex.what() << std::endl;
    return 1;
  }
  hmm_class.print(dir);

  if(search.empty()) {
    return 0;
  }

  fs::create_directories(tmp);

  if(!metacsst::subscan::run_scan_pipeline(search, out, tmp, options.threads, hmm_class)) {
    metacsst::app::cleanup_tmp_directory(tmp);
    return 1;
  }

  metacsst::app::cleanup_tmp_directory(tmp);
  return 0;
}
