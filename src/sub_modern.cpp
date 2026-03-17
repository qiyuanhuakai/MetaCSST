#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <thread>
#include <memory>
#include <filesystem>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <stdexcept>

///////////////////////////////////////////////////////////////////////////////////////////////
//Package: Metagenomic Complex Sequence Scanning Tool (MetaCSST)                             //
//Developer: Fazhe Yan                                                                       //
//Email: fazheyan33@163.com / ccwei@sjtu.edu                                                 //
//Department: Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University  //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "ghmm_modern.hpp"
using namespace std;
namespace fs = std::filesystem;

/*
This script is used to build a Weight Count Model according to the multiple sequence alignment result.In the meantime,get some conserved regions.And the build HMM models according to the motifs and using these HMM models to predict new structures in unknown sqeuences;
*/

/*
In the phase of finding motif,the training set sequences are clustered according to the result of multi-alignment.Foreach sub class,we find the best motif and build corresponding HMM model.When scaning for a new unknown sequence,these HMM models are combinded.This method may be not so efficient,but will be better in sensitivity as well as specificity.
*/

/*WorkFlow:
1>According to the trainging set,cluster the data to some sub clusters
2>Foreach cluster,find the best sequence motif using glam2 or muscle
3>Foreach motif,a GHMM model is built
4>All the GHMM models are used to scan for a new sequence
5>Combind the results of different GHMM model
*/

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
  else{
    int thread = 1;
    std::string search;
    std::string dir = "sbcsst_out";
    std::string config;
    for(int i = 0; i < argc; i++){
      if(strcmp(argv[i], "-thread") == 0 && i + 1 < argc)
        thread = atoi(argv[i + 1]);
      else if (strcmp(argv[i], "-in") == 0 && i + 1 < argc)
        search = argv[i + 1];
      else if (strcmp(argv[i], "-out") == 0 && i + 1 < argc)
        dir = argv[i + 1];
      else if (strcmp(argv[i], "-build") == 0 && i + 1 < argc)
        config = argv[i + 1];
      else if(strcmp(argv[i], "-h") == 0){
        metacsst::usage(argv[0]);
        return 0;
      }
    }
         
    if(config.empty()){
      metacsst::usage(argv[0]);
      return 0;
    }
    else{
      if(fs::exists(dir)){
        cout << "directory " << dir << " exists,it will be covered!" << endl;
        fs::remove_all(dir);
      }
      fs::create_directories(dir);
      
      std::string tmp = dir + "/tmp";
      std::string out = dir + "/out.txt";

      HMM_class hmm_class;
      try {
        hmm_class.init(config);
      } catch (const std::exception& ex) {
        cerr << "Config error: " << ex.what() << endl;
        return 1;
      }
      hmm_class.print(const_cast<char*>(dir.c_str()));

      if(!search.empty()){
        fs::create_directories(tmp);
        
        if(thread == 1){
          std::string out_tmp = tmp + "/out_tmp.txt";
          ScanArg arg;
          arg.search = search;
          arg.putout = out_tmp;
          arg.hmm_class = hmm_class;
          scanFile(arg);
          
          if(!fs::exists(out_tmp)){
            cerr << "Error: Temporary output file not created: " << out_tmp << endl;
            return 1;
          }
          
          if(out.empty()){
            std::ifstream infile(out_tmp);
            if(infile){
              cout << infile.rdbuf();
            }
          }
          else{
            fs::copy_file(out_tmp, out, fs::copy_options::overwrite_existing);
          }
        }
        else{
          int number = metacsst::split(search, thread, tmp);
          
          std::vector<std::thread> threads;
          std::vector<std::string> out_tmp(number);
          std::vector<std::string> sub_search(number);
          std::vector<std::unique_ptr<ScanArg>> args(number);
          
          for(int i = 0; i < number; i++){
            std::ostringstream oss_out, oss_search;
            oss_out << tmp << "/out_tmp_" << i << ".txt";
            out_tmp[i] = oss_out.str();
            
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
              cerr << "Error: Cannot create output file: " << out << endl;
              return 1;
            }
          }
          
          for(int i = 0; i < number; i++){
            std::ifstream infile(out_tmp[i]);
            if(infile){
              if(out.empty()){
                cout << infile.rdbuf();
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
      }
    }
  } 
  return 0;
}

void scanFile(ScanArg& arg){
  std::string name;
  name.reserve(100);
  
  std::ofstream out(arg.putout);
  if(!out){
    cerr << "Error: Cannot create output file: " << arg.putout << endl;
    return;
  }
  
  if(arg.search.empty()) return;
  
  std::ifstream in(arg.search);
  if(!in){
    cerr << "Error: Cannot open input file: " << arg.search << endl;
    return;
  }
  
  std::string line;
  while(std::getline(in, line)){
    metacsst::chomp(line);
    
    if(!line.empty() && line[0] == '>'){
      size_t space_pos = line.find(' ');
      size_t bracket_pos = line.find('[');
      size_t end_pos = std::min(space_pos, bracket_pos);
      if(end_pos == string::npos)
        name = line.substr(1);
      else
        name = line.substr(1, end_pos - 1);
    }
    else if(!line.empty()){
      vector<char> seq_buf(line.begin(), line.end());
      seq_buf.push_back('\0');
      
      struct OUT *result = arg.hmm_class.scanSeq(seq_buf.data());
      if(result->number > 0){
        out << name << "\n";
        out << line << "\n";
        for(int i = 0; i < result->number; i++){
          string seq_str(seq_buf.data());
          string matchSeq = metacsst::substr(seq_str, result->start[i], 
                                              result->end[i] - result->start[i] + 1);
          
          if(result->string[i] == 1)
            out << "Score:" << fixed << setprecision(2) << result->score[i] 
                << "\t+\tmatchSeq(" << result->start[i] << "-" << result->end[i] 
                << "):" << matchSeq << "\n";
          else{
            string matchSeq_complementary = metacsst::complementary(matchSeq);
            out << "Score:" << fixed << setprecision(2) << result->score[i] 
                << "\t-\tmatchSeq(" << result->start[i] << "-" << result->end[i] 
                << "):" << matchSeq_complementary << "\n";
          }
        }
      }
      free(result);
    }
  }
}
