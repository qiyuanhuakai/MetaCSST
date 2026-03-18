#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "ghmm_modern.hpp"

namespace fs = std::filesystem;

struct ScanArg {
  std::string search;
  std::string putout;
  HMM_class hmm_class;
};

void scanFile(ScanArg& argument);

int main(int argc,char* argv[]){
  if(argc < 3){
    metacsst::usage(argv[0]);
    return 0;
  }

  int thread = 1;
  std::string search;
  std::string dir = "sbcsst_out";
  std::string config;
  for(int i = 1; i < argc; ++i){
    const std::string_view opt(argv[i]);
    if(opt == "-thread" && i + 1 < argc) {
      try {
        thread = std::stoi(argv[++i]);
      } catch (const std::exception&) {
        metacsst::usage(argv[0]);
        return 1;
      }
    }
    else if (opt == "-in" && i + 1 < argc)
      search = argv[++i];
    else if (opt == "-out" && i + 1 < argc)
      dir = argv[++i];
    else if (opt == "-build" && i + 1 < argc)
      config = argv[++i];
    else if(opt == "-h"){
      metacsst::usage(argv[0]);
      return 0;
    }
  }

  if(thread <= 0){
    metacsst::usage(argv[0]);
    return 1;
  }

  if(config.empty()){
    metacsst::usage(argv[0]);
    return 0;
  }

  if(fs::exists(dir)){
    std::cout << "directory " << dir << " exists,it will be covered!" << std::endl;
    fs::remove_all(dir);
  }
  fs::create_directories(dir);

  const std::string tmp = dir + "/tmp";
  const std::string out = dir + "/out.txt";

  HMM_class hmm_class;
  try {
    hmm_class.init(config);
  } catch (const std::exception& ex) {
    std::cerr << "Config error: " << ex.what() << std::endl;
    return 1;
  }
  hmm_class.print(dir);

  if(search.empty()) {
    return 0;
  }

  fs::create_directories(tmp);

  if(thread == 1){
    const std::string out_tmp = tmp + "/out_tmp.txt";
    ScanArg arg;
    arg.search = search;
    arg.putout = out_tmp;
    arg.hmm_class = hmm_class;
    scanFile(arg);

    if(!fs::exists(out_tmp)){
      std::cerr << "Error: Temporary output file not created: " << out_tmp << std::endl;
      return 1;
    }

    if(out.empty()){
      std::ifstream infile(out_tmp);
      if(infile){
        std::cout << infile.rdbuf();
      }
    }
    else{
      fs::copy_file(out_tmp, out, fs::copy_options::overwrite_existing);
    }
  }
  else{
    const int number = metacsst::split(search, thread, tmp);

    std::vector<std::thread> threads;
    std::vector<std::string> out_tmp(number);
    std::vector<std::string> sub_search(number);
    std::vector<std::unique_ptr<ScanArg>> args(number);

    for(int i = 0; i < number; ++i){
      std::ostringstream oss_out;
      oss_out << tmp << "/out_tmp_" << i << ".txt";
      out_tmp[i] = oss_out.str();

      std::ostringstream oss_search;
      if(i >= 10)
        oss_search << tmp << "/split_" << i;
      else
        oss_search << tmp << "/split_0" << i;
      sub_search[i] = oss_search.str();

      args[i] = std::make_unique<ScanArg>();
      args[i]->search = sub_search[i];
      args[i]->putout = out_tmp[i];
      args[i]->hmm_class = hmm_class;

      threads.emplace_back(scanFile, std::ref(*args[i]));
    }

    for(auto& t : threads)
      t.join();

    std::ofstream outfile;
    if(!out.empty()){
      outfile.open(out);
      if(!outfile){
        std::cerr << "Error: Cannot create output file: " << out << std::endl;
        return 1;
      }
    }

    for(int i = 0; i < number; ++i){
      std::ifstream infile(out_tmp[i]);
      if(infile){
        if(out.empty()){
          std::cout << infile.rdbuf();
        }
        else{
          outfile << infile.rdbuf();
        }
      }
    }

    if(outfile.is_open())
      outfile.close();
  }
  fs::remove_all(tmp);
  return 0;
}

void scanFile(ScanArg& arg){
  std::string name;
  name.reserve(100);

  std::ofstream out(arg.putout);
  if(!out){
    std::cerr << "Error: Cannot create output file: " << arg.putout << std::endl;
    return;
  }

  if(arg.search.empty()) return;

  std::ifstream in(arg.search);
  if(!in){
    std::cerr << "Error: Cannot open input file: " << arg.search << std::endl;
    return;
  }

  std::string line;
  while(std::getline(in, line)){
    metacsst::chomp(line);

    if(!line.empty() && line[0] == '>'){
      const std::size_t space_pos = line.find(' ');
      const std::size_t bracket_pos = line.find('[');
      const std::size_t end_pos = std::min(space_pos, bracket_pos);
      if(end_pos == std::string::npos)
        name = line.substr(1);
      else
        name = line.substr(1, end_pos - 1);
    }
    else if(!line.empty()){
      auto result = arg.hmm_class.scanSeq(line);
      if(result->number > 0){
        out << name << "\n";
        out << line << "\n";
        const std::string& seq_str = line;
        for(int i = 0; i < result->number; i++){
          const std::string matchSeq = metacsst::substr(seq_str, result->start[i],
                                                       result->end[i] - result->start[i] + 1);

          if(result->string[i] == 1)
            out << "Score:" << std::fixed << std::setprecision(2) << result->score[i]
                << "\t+\tmatchSeq(" << result->start[i] << "-" << result->end[i]
                << "):" << matchSeq << "\n";
          else{
            const std::string matchSeq_complementary = metacsst::complementary(matchSeq);
            out << "Score:" << std::fixed << std::setprecision(2) << result->score[i]
                << "\t-\tmatchSeq(" << result->start[i] << "-" << result->end[i]
                << "):" << matchSeq_complementary << "\n";
          }
        }
      }
    }
  }
}
