#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <vector>
#include <string>
#include <filesystem>
#include <memory>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace fs = std::filesystem;

///////////////////////////////////////////////////////////////////////////////////////////////
//Package: Metagenomic Complex Sequence Scanning Tool (MetaCSST)                             //
//Developer: Fazhe Yan                                                                       //
//Email: fazheyan33@163.com / ccwei@sjtu.edu                                                 //
//Department: Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University  //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "ghmm_modern.hpp"
#include "fun_modern.hpp"
#include "config_modern.hpp"
using namespace std;

/*
This script is used to build a Weight Count Model according to the motif.In the meantime,get some conserved regions.And build GHMM models according to the motifs and using these GHMM models to predict new structures in unknown sqeuences;
*/

/*
In the phase of finding motif,the training set sequences are clustered according to the result of multi-alignment.Foreach sub class,we find the best motif and build corresponding GHMM model.When scaning for a new unknown sequence,these GHMM models are combinded.This method may be not so efficient,but will be better in sensitivity as well as specificity.
*/

/*WorkFlow:
1>According to the trainging set,cluster the data to some sub clusters
2>Foreach cluster,find the best sequence motif using glam2 or muscle
3>Foreach motif,a GHMM model is built
4>All the GHMM models are used to scan for a new sequence
5>Combind the results of different GHMM model
6>combind the results of different parts,for example:DGR=TR+VR+RT
*/

struct DGR {
  HMM_class TR; //clusters of TR
  HMM_class VR; //clusters of VR
  HMM_class RT; //clusters of RT
};

struct arg { 
  //arguments to scan the input file when using multi threads
  string search; //INPUT
  string putout; //OUTPUT
  SCAN dgrScan; 
  //the scan method,include three sub GHMM Models and a main GHMM model
};

void scanDGR(arg& argument); //scan the unknown sequence using the model
DGR buildDGRFromPath(const string& config_path);

int main(int argc,char* argv[]){
  if(argc < 3){
    usage(argv[0]);
    return 0;
  }
  else{
    int num_threads = 1; //thread number
    string config_path;
    string search; //unknown sequences file to scan
    string dir = "out_metacsst"; //out directory
    
    for(int i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-thread") == 0 && i + 1 < argc)
        num_threads = atoi(argv[i + 1]);
      else if(strcmp(argv[i], "-build") == 0 && i + 1 < argc)
        config_path = argv[i + 1];
      else if(strcmp(argv[i], "-in") == 0 && i + 1 < argc)
        search = argv[i + 1];
      else if(strcmp(argv[i], "-out") == 0 && i + 1 < argc)
        dir = argv[i + 1];
      else if(strcmp(argv[i], "-h") == 0){
        usage(argv[0]);
        return 0;
      }
    }
    
    if(!config_path.empty()){

      /*If the out directory exists,it will be covered*/
      if(fs::exists(dir)){
        cout << "directory " << dir << " exists,it will be covered!" << endl;
        fs::remove_all(dir);
      }
      
      /*mkdir: the out directory*/
      fs::create_directories(dir);
      
      /*tmp directory:used to save the temp results,including the split files and intermediate results*/
      string tmp = dir + "/tmp";
      /*out file:final result*/
      string out = dir + "/raw.gtf";

      DGR dgr;
      try {
        dgr = buildDGRFromPath(config_path);
      } catch (const std::exception& ex) {
        cerr << "Config error: " << ex.what() << endl;
        return 1;
      }
      //build DGR model accoring the config file
      SCAN dgrScan; //a parameter used in the multi-thread scaning
      dgrScan.init(dgr.TR, dgr.VR, dgr.RT, 10000);
      dgrScan.print(const_cast<char*>(dir.c_str()));

      if(!search.empty()){
        fs::create_directories(tmp);
        
        if(num_threads == 1){
          string out_tmp = tmp + "/out_tmp.txt";
          arg ARG;
          ARG.search = search;
          ARG.putout = out_tmp;
          ARG.dgrScan = dgrScan;

          scanDGR(ARG);
          
          // Concatenate output to final file
          if(!fs::exists(out)) {
            // Just copy the temp file to output location
            if(fs::exists(out_tmp)) {
              fs::copy(out_tmp, out, fs::copy_options::overwrite_existing);
            }
          } else {
            // Append temp file to existing output
            ifstream src(out_tmp);
            ofstream dst(out, ios::app);
            if(src && dst) {
              dst << src.rdbuf();
            }
          }
        }
        else{
          int number = split(search, num_threads, tmp); //split the input big file according to the number of threads
          vector<thread> threads; //threads
          vector<string> out_tmp(number);
          vector<string> sub_search(number);
          vector<arg> ARG(number); //arguments for each thread
          
          for(int i = 0; i < number; i++){
            ostringstream oss_out, oss_search;
            oss_out << tmp << "/out_tmp_" << i << ".txt";
            out_tmp[i] = oss_out.str();
            
            if(i >= 10)
              oss_search << tmp << "/split_" << i;
            else
              oss_search << tmp << "/split_" << setw(2) << setfill('0') << i;
            sub_search[i] = oss_search.str();
            
            ARG[i].search = sub_search[i];
            ARG[i].putout = out_tmp[i];
            ARG[i].dgrScan = dgrScan;
            
            threads.emplace_back(scanDGR, ref(ARG[i]));
          }

          // Wait for all threads to complete
          for(auto& t : threads){
            t.join();
          }
          
          // Concatenate all output files
          ofstream outfile(out);
          if(outfile) {
            for(int i = 0; i < number; i++){
              ifstream infile(out_tmp[i]);
              if(infile) {
                outfile << infile.rdbuf();
              }
            }
          }
        }
        
        // Clean up temp directory
        fs::remove_all(tmp);
      }
    }
  }
  return 0;
}


DGR buildDGRFromPath(const string& config_path){
  const auto dgr_cfg = metacsst::config::parse_dgr_motif_groups(config_path);
  DGR dgr;
  dgr.TR.init_groups(dgr_cfg.at("TR"));
  dgr.VR.init_groups(dgr_cfg.at("VR"));
  dgr.RT.init_groups(dgr_cfg.at("RT"));
  return dgr;
}

void scanDGR(arg& ARG){
  string tmp;
  ofstream out(ARG.putout);
  
  if(!ARG.search.empty() && out){
    ifstream IN(ARG.search);
    string name;
    
    while(getline(IN, tmp)){
      if(!tmp.empty() && tmp[0] == '>'){
        // Parse sequence name
        size_t end_pos = tmp.find_first_of(" \n");
        if(end_pos != string::npos && end_pos > 1) {
          name = tmp.substr(1, end_pos - 1);
        } else {
          name = tmp.substr(1);
        }
      }
      else if(!tmp.empty()){
        // Copy sequence to buffer for C-style scanning
        vector<char> seq_buf(tmp.begin(), tmp.end());
        seq_buf.push_back('\0');
        
        OUT* result = ARG.dgrScan.scanSeq(seq_buf.data());
        
        if(result->index == 1){ //index=1,there is sequence match
          for(int i = 0; i < result->number; i++){
            if(result->type[i] != 2){
              int start = result->start[i];
              int end = result->end[i];
              out << name << "\t";
              switch(result->type[i]){
                case 1: out << "TR\t"; break;
                case 3: out << "RT\t"; break;
              }
              
              string matchSeq = metacsst::substr(string(seq_buf.data()), start, end - start + 1);
              if(result->string[i] == 1)
                out << fixed << setprecision(2) << result->score[i] << "\t+\t" << start << "\t" << end << "\t" << matchSeq << "\n";
              else{
                string matchSeq_complementary = metacsst::complementary(matchSeq);
                out << fixed << setprecision(2) << result->score[i] << "\t-\t" << start << "\t" << end << "\t" << matchSeq_complementary << "\n";
              }
            }
          }
          string matchSeq = metacsst::substr(string(seq_buf.data()), result->total_start, result->total_end - result->total_start + 1);
          out << name << "\tDGR\t" << fixed << setprecision(2) << result->total_score << "\t*\t" << result->total_start << "\t" << result->total_end << "\t" << matchSeq << "\n";
        }
        free(result);
      }
    }
  }
}
