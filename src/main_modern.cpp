#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

#include "ghmm_modern.hpp"
#include "fun_modern.hpp"
#include "config_modern.hpp"

struct DGR {
  HMM_class TR;
  HMM_class VR;
  HMM_class RT;
};

struct ScanArgs {
  std::string search;
  std::string putout;
  SCAN dgrScan;
};

void scanDGR(ScanArgs& argument);
DGR buildDGRFromPath(const std::string& config_path);

int main(int argc,char* argv[]){
  if(argc < 3){
    metacsst::usage(argv[0]);
    return 0;
  }

  int num_threads = 1;
  std::string config_path;
  std::string search;
  std::string dir = "out_metacsst";

  for(int i = 1; i < argc; ++i) {
    const std::string_view opt(argv[i]);
    if(opt == "-thread" && i + 1 < argc) {
      try {
        num_threads = std::stoi(argv[++i]);
      } catch (const std::exception&) {
        metacsst::usage(argv[0]);
        return 1;
      }
    }
    else if(opt == "-build" && i + 1 < argc)
      config_path = argv[++i];
    else if(opt == "-in" && i + 1 < argc)
      search = argv[++i];
    else if(opt == "-out" && i + 1 < argc)
      dir = argv[++i];
    else if(opt == "-h"){
      metacsst::usage(argv[0]);
      return 0;
    }
  }

  if(num_threads <= 0) {
    metacsst::usage(argv[0]);
    return 1;
  }

  if(config_path.empty()) {
    metacsst::usage(argv[0]);
    return 0;
  }

  if(fs::exists(dir)){
    std::cout << "directory " << dir << " exists,it will be covered!" << std::endl;
    fs::remove_all(dir);
  }
  fs::create_directories(dir);

  const std::string tmp = dir + "/tmp";
  const std::string out = dir + "/raw.gtf";

  DGR dgr;
  try {
    dgr = buildDGRFromPath(config_path);
  } catch (const std::exception& ex) {
    std::cerr << "Config error: " << ex.what() << std::endl;
    return 1;
  }

  SCAN dgrScan;
  dgrScan.init(dgr.TR, dgr.VR, dgr.RT, 10000);
  dgrScan.print(dir);

  if(search.empty()) {
    return 0;
  }

  fs::create_directories(tmp);

  if(num_threads == 1){
    const std::string out_tmp = tmp + "/out_tmp.txt";
    ScanArgs args;
    args.search = search;
    args.putout = out_tmp;
    args.dgrScan = dgrScan;

    scanDGR(args);

    if(!fs::exists(out)) {
      if(fs::exists(out_tmp)) {
        fs::copy(out_tmp, out, fs::copy_options::overwrite_existing);
      }
    } else {
      std::ifstream src(out_tmp);
      std::ofstream dst(out, std::ios::app);
      if(src && dst) {
        dst << src.rdbuf();
      }
    }
  }
  else{
    const int number = metacsst::split(search, num_threads, tmp);
    std::vector<std::thread> threads;
    std::vector<std::string> out_tmp(number);
    std::vector<std::string> sub_search(number);
    std::vector<ScanArgs> args(number);

    for(int i = 0; i < number; ++i){
      std::ostringstream oss_out;
      oss_out << tmp << "/out_tmp_" << i << ".txt";
      out_tmp[i] = oss_out.str();

      std::ostringstream oss_search;
      if(i >= 10)
        oss_search << tmp << "/split_" << i;
      else
        oss_search << tmp << "/split_" << std::setw(2) << std::setfill('0') << i;
      sub_search[i] = oss_search.str();

      args[i].search = sub_search[i];
      args[i].putout = out_tmp[i];
      args[i].dgrScan = dgrScan;

      threads.emplace_back(scanDGR, std::ref(args[i]));
    }

    for(auto& t : threads){
      t.join();
    }

    std::ofstream outfile(out);
    if(outfile) {
      for(int i = 0; i < number; ++i){
        std::ifstream infile(out_tmp[i]);
        if(infile) {
          outfile << infile.rdbuf();
        }
      }
    }
  }

  fs::remove_all(tmp);
  return 0;
}

DGR buildDGRFromPath(const std::string& config_path){
  const auto dgr_cfg = metacsst::config::parse_dgr_motif_groups(config_path);
  DGR dgr;
  dgr.TR.init_groups(dgr_cfg.at("TR"));
  dgr.VR.init_groups(dgr_cfg.at("VR"));
  dgr.RT.init_groups(dgr_cfg.at("RT"));
  return dgr;
}

void scanDGR(ScanArgs& args){
  std::string tmp_line;
  std::ofstream out(args.putout);

  if(!args.search.empty() && out){
    metacsst::LineReader in(args.search);
    if(!in.is_open()) {
      std::cerr << "Error: Cannot open input file: " << args.search << std::endl;
      return;
    }
    std::string name;

    while(in.getline(tmp_line)){
      if(!tmp_line.empty() && tmp_line[0] == '>'){
        const std::size_t end_pos = tmp_line.find_first_of(" \n");
        if(end_pos != std::string::npos && end_pos > 1) {
          name = tmp_line.substr(1, end_pos - 1);
        } else {
          name = tmp_line.substr(1);
        }
      }
      else if(!tmp_line.empty()){
        auto result = args.dgrScan.scanSeq(tmp_line);

        if(result->index == 1){
          const std::string& seq_str = tmp_line;
          for(int i = 0; i < result->number; ++i){
            if(result->type[i] != 2){
              const int start = result->start[i];
              const int end = result->end[i];
              out << name << "\t";
              switch(result->type[i]){
                case 1: out << "TR\t"; break;
                case 3: out << "RT\t"; break;
              }

              const std::string match_seq = metacsst::substr(seq_str, start, end - start + 1);
              if(result->string[i] == 1)
                out << std::fixed << std::setprecision(2) << result->score[i] << "\t+\t" << start << "\t" << end << "\t" << match_seq << "\n";
              else{
                const std::string match_seq_complementary = metacsst::complementary(match_seq);
                out << std::fixed << std::setprecision(2) << result->score[i] << "\t-\t" << start << "\t" << end << "\t" << match_seq_complementary << "\n";
              }
            }
          }
          const std::string match_seq = metacsst::substr(seq_str, result->total_start, result->total_end - result->total_start + 1);
          out << name << "\tDGR\t" << std::fixed << std::setprecision(2) << result->total_score << "\t*\t" << result->total_start << "\t" << result->total_end << "\t" << match_seq << "\n";
        }
      }
    }
  }
}
