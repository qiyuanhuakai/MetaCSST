#ifndef GHMM_MODERN_HPP
#define GHMM_MODERN_HPP

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <array>
#include <string>
#include <string_view>
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <memory>

///////////////////////////////////////////////////////////////////////////////////////////////
//Package: Metagenomic Complex Sequence Scanning Tool (MetaCSST)                             //
//Developer: Fazhe Yan                                                                       //
//Email: fazheyan33@163.com / ccwei@sjtu.edu                                                 //
//Department: Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University  //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "fun_modern.hpp"
#include "config_modern.hpp"
using std::array;
using std::exception;
using std::ifstream;
using std::string;
using std::to_string;
using std::vector;

//=============================================================================================
// Data structures for GHMM (Generalized Hidden Markov Model)
//=============================================================================================

/**
 * @brief Pattern box derived from multiple sequence alignment.
 *
 * Represents a conserved motif region. Stores the position-specific
 * scoring matrix (PSSM) and alignment counts for A/T/C/G.
 */
struct pattern {
  int length; //pattern box length
  vector<int> sum; //the sum of A,T,C,G appearance time in every position
  vector<vector<int>> matrix; //alignment metrix
  vector<vector<float>> score; //Position specific scoring matrix
  float max; //max score of a random sequence
  float min; //min score of a random sequence
  int pos_start; //start position of this pattern box in the training set
  int pos_end; //end position pf the pattern box in the training set
};

/**
 * @brief HMM state box.
 *
 * Each box corresponds to one state in the GHMM and is represented by
 * a PSSM derived from the training alignment.
 */
struct box { //every box is a state in the Hidden Markov Model
  int length; //state matrix length
  int pos_start; //start position of the pattern box in the training set sequence
  int pos_end; //end position of the pattern box in the training set sequence
  float max,min; //max score and min score
  vector<vector<int>> align; //align matrix of the state box
  vector<vector<float>> score; //scoring matrix of the state box
};

/**
 * @brief Single match state returned by a GHMM scan.
 */
struct match_state {
  int start;
  int end;
  float score;
  int strand;
};

/**
 * @brief Sub-structure match identified by a single HMM.
 *
 * Used for TR (Template Repeat), VR (Variable Repeat) or RT
 * (Reverse Transcriptase) motifs.
 */
struct sub_hmm { //sub HMM strctures ,such as TR/VR/RT
  int start; //start site of the sub HMM in the input sequence
  int end; //end site
  float score; //match score of this sub HMM structure
  int index; //index=-1 -> init;  index=0->TR;  index=1->VR  index=2->RT;
};

/**
 * @brief Container for scanning results.
 *
 * Holds either sub-HMM matches (TR/VR/RT) or a complete DGR
 * (Diversity-Generating Retroelement) detection result.
 */
struct out_result { //scaning result,for sub HMM(TR/VR/RT) or the total DGR
  int number; //matchSeq number
  std::array<float, S> score{}; //score of the matches
  std::array<int, S> start{}; //start site
  std::array<int, S> end{}; //end site
  std::array<int, S> string{}; //string of the match,1::'+' or 2::'-'

  std::array<int, S> type{}; //used only for scaning fot total DGR;1->TR,2->VR,3->RT;
  int total_start; //used only for scaning for DGR
  int total_end; //used only for scaning for DGR
  float total_score; //used only for scaning for DGR
  int index; //used only for DGR;0->no hot,1->DGR in found
};

//=============================================================================================
// Single GHMM model (hmm_model)
//   - Initializes from training alignment
//   - Scans sequences using the Viterbi algorithm on PSSM states
//=============================================================================================

class hmm_model{
 private:
  int len; //length cuttof used to ensure a state
  int size; //state number in the Hidden Markov Model
  int window; //max box length,used as a window size
  int gap; //max gap length between two states
  float cutoff; //cutoff value for a subsequence matching a state box(0~1)
  vector<string> name; //state names
  vector<box> state; //every state is a box,represented by a scoring metrix(Position Specific Scoring Metrix)
  vector<vector<float>> trans; //Transition probability Metrix beteen the states
  vector<float> start,end; //start and end probability for the states
 public:
  float score_cutoff; //score cutoff for a sequence belong to this HMM model

  hmm_model(): len(0), size(0), window(0), gap(0), cutoff(0.0f), score_cutoff(0.0f) {}

  hmm_model(const vector<vector<float>>& transition,
      const vector<pattern>& metrix,
      int number,
      float value,
      int window_size,
      int gap_length,
      int state_length,
      float seq_score_cutoff){
    init(transition,metrix,number,value,window_size,gap_length,state_length,seq_score_cutoff);
  }

  void init(const vector<vector<float>>& transition,
            const vector<pattern>& metrix,
            int number,
            float cutoff_value,
            int window_size,
            int gap_length,
            int state_length,
            float seq_score_cutoff);
  void print(const std::string& dir);
  std::unique_ptr<out_result> scan_seq_single(std::string_view seq); //only scan for positive string
  std::unique_ptr<out_result> scan_seq_full(std::string_view seq); //scan for the both two directions

  std::unique_ptr<out_result> scanSeqSingle(std::string_view seq) { return scan_seq_single(seq); }
  std::unique_ptr<out_result> scanSeqFull(std::string_view seq) { return scan_seq_full(seq); }
};

inline bool arg_equals(std::string_view arg, std::string_view option) {
  return arg == option;
}

void hmm_model::init(const vector<vector<float>>& transition,
               const vector<pattern>& metrix,
               int number,
               float cutoff_value,
               int window_size,
               int gap_length,
               int state_length,
               float seq_score_cutoff){
  len = state_length;
  gap = gap_length;
  size = number;
  window = window_size;
  cutoff = cutoff_value;
  score_cutoff = seq_score_cutoff;

  name.assign(size,string());
  state.assign(size,box());
  trans.assign(size,vector<float>(size,0.0f));
  start.assign(size,0.0f);
  end.assign(size,0.0f);

  for(int i=0;i<size;i++)
    name[i] = "pattern" + to_string(i);

  for(int i=0;i<size+2;i++) //get the Transition probability Metrix
    for(int j=0;j<size+2;j++)
      if(i==0 && j>=1 && j<=size)
        start[j-1] = transition[i][j];
      else if(j==size+1 && i>=1 && i<=size)
        end[i-1] = transition[i][j];
      else if(i>=1 && i<=size && j>=1 && j<=size)
        trans[i-1][j-1]=transition[i][j];

  for(int i=0;i<size;i++){ //get the states accoring the input patterns
    state[i].length = metrix[i].length;
    state[i].pos_start = metrix[i].pos_start;
    state[i].pos_end = metrix[i].pos_end;
    state[i].max = metrix[i].max;
    state[i].min = metrix[i].min;

    state[i].align.assign(4,vector<int>(state[i].length,0));
    state[i].score.assign(4,vector<float>(state[i].length,0.0f));
    for(int j=0;j<4;j++)
      for(int k=0;k<metrix[i].length;k++){
        state[i].align[j][k] = metrix[i].matrix[j][k];
        state[i].score[j][k] = metrix[i].score[j][k];
      }
  }
}

/**
 * @brief Persist the trained model matrices to disk.
 *
 * Writes alignment matrices to "align.txt" and scoring / transition
 * probabilities to "score.txt" under the given directory.
 */
void hmm_model::print(const std::string& dir){

  const auto align_path = std::filesystem::path(dir) / "align.txt";
  const auto score_path = std::filesystem::path(dir) / "score.txt";

  std::ofstream fp1(align_path, std::ios::app);
  std::ofstream fp2(score_path, std::ios::app);
  if(!fp1 || !fp2){
    return;
  }

  for(int i=0;i<size;i++)
    if(state[i].length >= len){
      fp1 << "align matrix:\n";
      for(int j=0;j<4;j++){
        switch(j){
        case 0:fp1 << "A\t";break;
        case 1:fp1 << "T\t";break;
        case 2:fp1 << "C\t";break;
        case 3:fp1 << "G\t";break;
        }

        for(int k=0;k<state[i].length;k++)
          if(k == state[i].length -1)
            fp1 << state[i].align[j][k] << '\n';
          else
            fp1 << state[i].align[j][k] << '\t';
      }
      fp2 << "scoring matrix:\n";
      for(int j=0;j<4;j++){
        switch(j){
        case 0:fp2 << "A\t";break;
        case 1:fp2 << "T\t";break;
        case 2:fp2 << "C\t";break;
        case 3:fp2 << "G\t";break;
        }

        for(int k=0;k<state[i].length;k++)
          if(k == state[i].length -1)
            fp2 << std::fixed << std::setprecision(2) << state[i].score[j][k] << '\n';
          else
            fp2 << std::fixed << std::setprecision(2) << state[i].score[j][k] << '\t';
      }
    }

  fp2 << "Start Probabilities\t";
  for(int i=0;i<size;i++)
    if(state[i].length >= len)
      fp2 << std::fixed << std::setprecision(2) << start[i] << '\t';
  fp2 << "\nEnding Probabilities\t";
  for(int i=0;i<size;i++)
    if(state[i].length >= len)
      fp2 << std::fixed << std::setprecision(2) << end[i] << '\t';
  fp2 << "\nTrasnition Probability Matrix\n";
  for(int i=0;i<size;i++)
    if(state[i].length >= len){
      for(int j=0;j<size;j++)
        if(state[j].length >= len)
          fp2 << std::fixed << std::setprecision(2) << trans[i][j] << '\t';
      fp2 << '\n';
    }

  fp2 << "Window size:" << window << '\n';
  fp2 << "Gap length:" << gap << '\n';
  fp2 << "Match score cutoff:" << std::fixed << std::setprecision(2) << score_cutoff << '\n';
  fp1 << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  fp2 << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}


/**
 * @brief Scan a single DNA strand for motif matches using the Viterbi algorithm.
 *
 * Workflow:
 *   1) Slide a window across the sequence and score each sub-sequence
 *      against every state PSSM. Record high-scoring motif hits.
 *   2) Split the hit list into independent search spaces whenever the
 *      distance between consecutive hits exceeds @c gap.
 *   3) Inside each search space, apply the Viterbi algorithm to find the
 *      maximum-likelihood path through the GHMM states.
 *   4) Emit all paths whose total score exceeds @c score_cutoff.
 */
std::unique_ptr<out_result> hmm_model::scan_seq_single(std::string_view seq){
  /*workflow of scaning for TR,VR or RT:
    (1)foreach sub state,set down the matching subSeqs(position as well as score).
    (2)according to the gap length,the whole sequence will be splitted to some search space.
    (3)foreach search space,based on the subSeqs,build many paths according to the veterbi algorithm.
    Every path is a solution for the problem,and the path with the best score will be choosed.
    (4)save all the result satisfying the requirements:score(path) > score_cutoff
   */

  auto result = std::make_unique<out_result>();
  result->number=0;

  const int S2 = static_cast<int>(seq.size());

  vector<int> state_pos(S2,-1);
  vector<int> state_index(S2,-1);
  vector<float> state_score(S2,0.0f);
  int number = 0;

  /*Motif scaning strategy:
    if there is a sequence matched,then skip length will be the length of this motif rather then 1.This methos is able to accelerate the process,but may miss the best match in the mean while,which reduces the whole score,leading to a wrong result.
  */

  //-------------------------------------------------------------------------------------------
  // Stage 1: Motif scanning - score every window against each state's PSSM.
  //-------------------------------------------------------------------------------------------
  if(S2 > window){
    for(int i=0;i<S2-window-1;){ //scan the sequence for every stat,set down the position as well as matching score
      int index=0;
      //index:in this position,exists a state match? 0:no  1:yes
      for(int j=0;j<size;j++){
        if(state[j].length >= len) {
          float tmp=0.0f;

          for(int k=0;k < state[j].length;k++){
            switch(seq[static_cast<std::size_t>(i+k)]){
            case 'A':tmp += state[j].score[0][k];break;
            case 'T':tmp += state[j].score[1][k];break;
            case 'C':tmp += state[j].score[2][k];break;
            case 'G':tmp += state[j].score[3][k];break;
            case 'a':tmp += state[j].score[0][k];break;
            case 't':tmp += state[j].score[1][k];break;
            case 'c':tmp += state[j].score[2][k];break;
            case 'g':tmp += state[j].score[3][k];break;
            default:tmp += state[j].min/state[j].length;
              //score for base 'N'
            }
          }

          if(tmp > cutoff*state[j].max && tmp > 0){
            state_pos[number] = i;
            state_index[number] = j;
            state_score[number] = tmp;
            i += state[j].length;
            number++;
            index = 1;
            break;
          }
        }
      }
      if(index == 0)
        i++;
    }
  }

  /*There may be many mstchSeqs for an input sequence,if gap length more than the gap,then searching for a new matchSeq*/

  vector<int> search_start(S2,-1);
  vector<int> search_end(S2,-1);

  //-------------------------------------------------------------------------------------------
  // Stage 2: Divide hits into independent search spaces using the gap threshold.
  //-------------------------------------------------------------------------------------------
  int search_number=0; //according to the gap length,devide the total search space to some small search space
  for(int i=0;i<number;i++){
    if(i==0)
      search_start[0] = 0;
    if(i<number-1 && state_pos[i+1]-state_pos[i] > gap){ //the gap length is more than gap
      search_end[search_number] = i;
      search_number++;
      search_start[search_number] = i+1;
    }
    if(i==number-1){
      search_end[search_number] = number-1;
      search_number++;
    }
  }

  /*For each search space,using veterbi algorithm to find the best path.if the path score is larger than the scpre cuttof,save the score and start position as well as the end position*/

  vector<float> veterbi_score(number,0.0f); //score of the possible path using veterbi algorithm
  vector<int> veterbi_start(number,-1); //matching start of the veterbi path

  //-------------------------------------------------------------------------------------------
  // Stage 3 & 4: Viterbi dynamic programming per search space + emit results.
  //-------------------------------------------------------------------------------------------
  int num=0; //number of matchSeqs in all the search space
  for(int k=0;k<search_number;k++){
    for(int i=search_start[k];i <= search_end[k];i++){
      float max = 0.0f;

      if(start[state_index[i]] != 0){
        float score_tmp = start[state_index[i]] * state_score[i]; //start probability
        if(score_tmp > max){
          max = score_tmp;
          veterbi_start[i] = state_pos[i];
        }
      }
      if(i > search_start[k])
        for(int j=search_start[k];j<i;j++){
          float score_tmp = 0.0f;
          if(trans[state_index[j]][state_index[i]] != 0 && veterbi_score[j] != 0){
            score_tmp = veterbi_score[j] + state_score[i] * trans[state_index[j]][state_index[i]];  //state transition
            if(score_tmp > max){
              max = score_tmp;
              veterbi_start[i] = veterbi_start[j];
            }
          }
        }
      veterbi_score[i] = max;
    }

    float score_search=0.0f;int start_search=-1;int end_search=-1;
    /*for each search space,there may be many pathway if exists more then one state.According to the veterbi algorithm,choose the path with the highest score*/
    for(int j=search_start[k];j<=search_end[k];j++){
      if(end[state_index[j]] != 0 && veterbi_score[j] != 0){
        float score_tmp = veterbi_score[j] + end[state_index[j]];
        if(score_tmp > score_search){
          score_search = score_tmp;
          end_search = state_pos[j]+state[state_index[j]].length-1;
          start_search = veterbi_start[j];
        }
      }
    }
    if(score_search > score_cutoff){
      result->score[num] = score_search;
      result->start[num] = start_search;
      result->end[num] = end_search;
      num++;
      result->number ++;
    }
  }

  return result;
}

/**
 * @brief Scan both DNA strands (positive + reverse complement).
 *
 * Coordinates on the reverse strand are converted back to the
 * original (positive) orientation before merging.
 */
std::unique_ptr<out_result> hmm_model::scan_seq_full(std::string_view seq){
  auto result = std::make_unique<out_result>();
  result->number=0;

  auto result1 = scan_seq_single(seq); //positive chain

  string seq_complementary = metacsst::complementary(string(seq));
  auto result2 = scan_seq_single(seq_complementary);

  const int length = static_cast<int>(seq.size());
  int num=0;
  if(result1->number > 0)
    for(int i=0;i<result1->number;i++,num++){
      result->score[num] = result1->score[i];
      result->start[num] = result1->start[i];
      result->end[num] = result1->end[i];
      result->string[num] = 1;
    }

  if(result2->number > 0)
    for(int i=0;i<result2->number;i++,num++){
      result->score[num] = result2->score[i];
      result->start[num] = length - result2->end[i];
      result->end[num] = length - result2->start[i];
      result->string[num] = 2;
    }

  result->number = num;
  return result;
}

//=============================================================================================
// HMM builder from multiple alignment file
//=============================================================================================

/**
 * @brief Build a single GHMM from a training alignment file.
 *
 * Steps:
 *   1) Parse CLI arguments (coverage, length, score ratio, gap, etc.).
 *   2) Read the aligned sequences and compute per-position counts of A/T/C/G.
 *   3) Derive the Position-Specific Scoring Matrix (PSSM) and Information Content.
 *   4) Extract conserved pattern boxes using the coverage cutoff.
 *   5) Count state transitions in the training set and build the transition matrix.
 *   6) Compute an empirical score cutoff from the training sequences.
 */
hmm_model build_hmm(const std::vector<std::string>& args){
  int L=0; //sequence length in the multiAlignment result
  float cov=0.9f; //coverage cuttof in every position to make a pattern box
  int box_len_cutoff=7; //pattern box length cutoff
  int max_box_length=0; //max length of the insured patterns
  int pattern_number=0; //pattern box number in total
  float state_score=0.4f; //cuttof of the score(ratio) used to ensure a state
  float ratio=1.0f; //TP value to control when testing
  int gap = 100; //gap between the states
  float ic = 0.5f; //IC value cuttof

  hmm_model hmm;

  if(args.size() < 2){
    usage(args.empty() ? std::string("hmm") : args.front());
  }
  else{
    int i,j,k;
    string in_path;
    for(i=0;i<static_cast<int>(args.size());i++)
      if(arg_equals(args[static_cast<std::size_t>(i)],"-build") && i+1<static_cast<int>(args.size()))
        in_path = args[static_cast<std::size_t>(i+1)];
      else if(arg_equals(args[static_cast<std::size_t>(i)],"-cov") && i+1<static_cast<int>(args.size()))
        cov = std::stof(args[static_cast<std::size_t>(i+1)]);
      else if (arg_equals(args[static_cast<std::size_t>(i)],"-len") && i+1<static_cast<int>(args.size()))
        box_len_cutoff = std::stoi(args[static_cast<std::size_t>(i+1)]);
      else if (arg_equals(args[static_cast<std::size_t>(i)],"-score") && i+1<static_cast<int>(args.size()))
        state_score = std::stof(args[static_cast<std::size_t>(i+1)]);
      else if (arg_equals(args[static_cast<std::size_t>(i)],"-ratio") && i+1<static_cast<int>(args.size()))
        ratio = std::stof(args[static_cast<std::size_t>(i+1)]);
      else if (arg_equals(args[static_cast<std::size_t>(i)],"-gap") && i+1<static_cast<int>(args.size()))
        gap = std::stoi(args[static_cast<std::size_t>(i+1)]);
      else if (arg_equals(args[static_cast<std::size_t>(i)],"-ic") && i+1<static_cast<int>(args.size()))
        ic = std::stof(args[static_cast<std::size_t>(i+1)]);
      else if(arg_equals(args[static_cast<std::size_t>(i)],"-h"))
        usage(args.front());

    if(in_path.empty() || cov <= 0 || cov >1){
      usage(args.front());
      return hmm;
    }
    else{
      (void)ic;

      ifstream in_first(in_path);
      if(!in_first.is_open()){
        usage(args.front());
        return hmm;
      }

      string tmp;
      while(getline(in_first,tmp)){
        if(!tmp.empty() && tmp.back() == '\r')
          tmp.pop_back();
        if(judge(tmp) == 0)
          if(L == 0)
            L = static_cast<int>(tmp.length()); //get the aligned sequence length
      }
      in_first.close();

      if(L <= 0)
        return hmm;

      ifstream in(in_path);
      if(!in.is_open())
        return hmm;

      //---------------------------------------------------------------------------------------
      // Build the alignment count matrix and background prior probabilities.
      //---------------------------------------------------------------------------------------
      vector<vector<int>> count(4,vector<int>(L,0)); //align metrix,store the appearance times of A,T,C,G in every position
      vector<float> priori(4,0.0f); //priori probability of the background
      vector<int> pri_tmp(4,0); //a buffer to store the appearance times of A,T,C,G,used to calculate the  priori probability
      vector<int> sum(L,0); //sum of the symbol appearance in every position
      vector<float> IC(L,0.0f); //information value in every position
      int number=0; //sequence number

      while(getline(in,tmp)){  //read the input sequences and build the alignment matrix
        if(!tmp.empty() && tmp.back() == '\r')
          tmp.pop_back();
        if(judge(tmp) == 0){
          number++;
          int line_len = static_cast<int>(tmp.length());
          for(i=0;i<line_len && i<L;i++)
            switch(tmp[i]){
            case 'A':count[0][i]++;break;
            case 'T':count[1][i]++;break;
            case 'C':count[2][i]++;break;
            case 'G':count[3][i]++;break;
            }
        }
      }

      for(i=0;i<L;i++){
        for(j=0;j<4;j++)
          pri_tmp[j] += count[j][i];
        sum[i]=count[0][i]+count[1][i]+count[2][i]+count[3][i];
      }

      int sum_pri = pri_tmp[0]+pri_tmp[1]+pri_tmp[2]+pri_tmp[3];
      if(sum_pri == 0 || number == 0){
        in.close();
        return hmm;
      }

      for(j=0;j<4;j++)
        priori[j] = pri_tmp[j]*1.0f/sum_pri;

      //---------------------------------------------------------------------------------------
      // Compute Position-Specific Scoring Matrix (PSSM) and Information Content.
      //---------------------------------------------------------------------------------------
      vector<vector<float>> score(5,vector<float>(L,0.0f)); //Positon Specific Scoring Metrix,including the positon coverage(score[0])
      for(i=0;i<L;i++){
        score[0][i]=sum[i]*1.0f/number; //calculate position coverage
        IC[i] = 0.0f;

        for(j=0;j<4;j++)
          if(count[j][i] != 0)
            IC[i] += static_cast<float>(log(count[j][i]/(priori[j]*sum[i]))*count[j][i]/sum[i]); //calculate IC value
      }
      for(i=0;i<4;i++)
        for(j=0;j<L;j++)
          if(count[i][j] == 0)
            score[i+1][j]=-10.0f;
          else
            score[i+1][j] = static_cast<float>(log(count[i][j]/(priori[i]*sum[j]))*sum[j]/number); //formula to calculate the score of every position

      //---------------------------------------------------------------------------------------
      // Extract conserved pattern boxes based on the coverage threshold.
      //---------------------------------------------------------------------------------------
      vector<pattern> scan(M);  //store some subPattern in the whole scoring matrix
      for(std::size_t scan_idx = 0; scan_idx < scan.size(); ++scan_idx){
        scan[scan_idx].length=0;
        scan[scan_idx].max=0.0f;
        scan[scan_idx].min=0.0f;
        scan[scan_idx].pos_start=-1;
        scan[scan_idx].pos_end=-1;
        scan[scan_idx].score.assign(4,vector<float>(P,0.0f));
        scan[scan_idx].matrix.assign(4,vector<int>(P,0));
        scan[scan_idx].sum.assign(P,0);
      }

      for(i=0,j=0,k=0;i<L;i++){
        if(score[0][i] >= cov){
          //  if(score[0][i] >= cov && IC[i] >= ic){
          if(scan[k].pos_start == -1)
            scan[k].pos_start = i;
          scan[k].length += 1;
          scan[k].score[0][j]=score[1][i];scan[k].matrix[0][j]=count[0][i];
          scan[k].score[1][j]=score[2][i];scan[k].matrix[1][j]=count[1][i];
          scan[k].score[2][j]=score[3][i];scan[k].matrix[2][j]=count[2][i];
          scan[k].score[3][j]=score[4][i];scan[k].matrix[3][j]=count[3][i];
          scan[k].sum[j]=sum[i];
          j++;
          if(i==L-1){
            scan[k].pos_end=i;
            k++;
          }
        }
        else if(i>0 && score[0][i-1] >= cov){
          //else if(i>0 && score[0][i-1] >= cov && IC[i-1] >= ic){
          scan[k].pos_end=i-1;
          k++;
          j=0;
        }
      }

      pattern_number = k; //pattern  number satisfying given conditions(such as coverage,length,i.e.)
      //calculate the max_box_length,min/max score,Information content and corresponding p value of every pattern
      for(i=0;i<pattern_number;i++){
        if(scan[i].length > max_box_length)
          max_box_length = scan[i].length;
        for(j=0;j < scan[i].length;j++){
          float max=scan[i].score[0][j];
          float min=scan[i].score[0][j];
          for(k=1;k<=3;k++){
            if(scan[i].score[k][j] > max)
              max = scan[i].score[k][j];
            if(scan[i].score[k][j] < min)
              min = scan[i].score[k][j];
          }
          scan[i].max += max;
          scan[i].min += min;
        }
      }

      //---------------------------------------------------------------------------------------
      // Build the state transition probability matrix from the training sequences.
      //---------------------------------------------------------------------------------------
/* buildTransition */
      int num=pattern_number+2; //state number:pattern_number+2(start+patterns+end)
      vector<vector<int>> trans_count(num,vector<int>(num,0)); //state transition frequency(times) in the training set
      vector<vector<float>> transition(num,vector<float>(num,0.0f)); //start transition probability metrix

      in.close();
      in.open(in_path);
      if(!in.is_open())
        return hmm;

      while(getline(in,tmp)){
        if(!tmp.empty() && tmp.back() == '\r')
          tmp.pop_back();
        if(judge(tmp) != -1){
          vector<int> pattern_index(pattern_number,-1);
          for(i=0,k=0;i<pattern_number;i++){
            if(scan[i].length >= box_len_cutoff){
              float score_tmp=0.0f;
              string sub = metacsst::substr(tmp,scan[i].pos_start,scan[i].length);
              for(j=0;j<static_cast<int>(sub.length());j++)
                switch(sub[j]){
                case 'A':score_tmp += scan[i].score[0][j];break;
                case 'T':score_tmp += scan[i].score[1][j];break;
                case 'C':score_tmp += scan[i].score[2][j];break;
                case 'G':score_tmp += scan[i].score[3][j];break;
                default:score_tmp += 0;
                }
              if(score_tmp > state_score*scan[i].max){
                pattern_index[k] = i;
                k++;
              }
            }
          }
          for(i=0;i<k;i++){
            int id = pattern_index[i];
            if(i == 0)
              trans_count[0][id+1]++;
            if(i == k-1)
              trans_count[id+1][num-1]++;
            if(i>0){
              int id2 = pattern_index[i-1];
              trans_count[id2+1][id+1]++;
            }
          }
        }
      }

      for(i=0;i<num;i++){
        int sum_line = 0;
        for(j=0;j<num;j++)
          sum_line += trans_count[i][j];
        for(j=0;j<num;j++)
          transition[i][j] = sum_line>0?trans_count[i][j]*1.0f/sum_line:0.0f; //calculate probability
      }
      in.close();
/* BuildTransition End */

      //---------------------------------------------------------------------------------------
      // Derive a score cutoff from the training set match scores.
      //---------------------------------------------------------------------------------------
/*Get the score cottof for a HMM matching according to the training set */
      int line = 0;
      vector<float> score_train;
      score_train.reserve(1024);

      in.open(in_path); //get the cuttof value of scaning according to the annotated training set
      if(!in.is_open())
        return hmm;

      while(getline(in,tmp)){
        if(!tmp.empty() && tmp.back() == '\r')
          tmp.pop_back();
        if(judge(tmp) == 0){
          int pri = -1; //the previous state
          float path_score=0.0f;
          for(i=0;i<pattern_number;i++){
            if(scan[i].length >= box_len_cutoff){
              float score_tmp=0.0f;
              string sub = metacsst::substr(tmp,scan[i].pos_start,scan[i].length);
              for(j=0;j<static_cast<int>(sub.length());j++)
                switch(sub[j]){
                case 'A':score_tmp += scan[i].score[0][j];break;
                case 'T':score_tmp += scan[i].score[1][j];break;
                case 'C':score_tmp += scan[i].score[2][j];break;
                case 'G':score_tmp += scan[i].score[3][j];break;
                default:score_tmp += 0;
                }
              if(score_tmp > state_score*scan[i].max){
                if(pri == -1)
                  path_score += score_tmp*transition[0][i+1];
                else
                  path_score += score_tmp*transition[pri+1][i+1];
                pri = i;
              }
            }
          }
          if(pri != -1)
            path_score += transition[pri+1][pattern_number+1];
          score_train.push_back(path_score);
          line++;
        }
      }
      in.close();

      score_train.push_back(0.0f);
      float *score_train_ptr = score_train.data();
      float score_cutoff = cutoff(&score_train_ptr, line, ratio);

      hmm = hmm_model(transition,scan,pattern_number,state_score,max_box_length,gap,box_len_cutoff,score_cutoff);
    }
  }
  return hmm;
}

//=============================================================================================
// HMM class (hmm_class)
//   - Manages a cluster of hmm_model instances belonging to the same motif class
//=============================================================================================

class hmm_class { //clusters of GHMM model
 public:
  vector<hmm_model> hmm;
  //multi similar GHMM models,which belongs to different classes
  int _number; //number of clusters(models)

  hmm_class(): _number(0) {}
  explicit hmm_class(const string& config): _number(0) { init(config); }

  void init(const string& config); //initialization,based on the config file
  void init_groups(const std::vector<metacsst::config::OrderedKeyValues>& motif_groups);
  void print(const std::string& dir);
  std::unique_ptr<out_result> scan_seq(std::string_view seq); //scaning a new sequence

  std::unique_ptr<out_result> scanSeq(std::string_view seq) { return scan_seq(seq); }
};

/**
 * @brief Initialize from a list of motif groups.
 *
 * Each group supplies command-line-style arguments that are forwarded
 * to build_hmm() to construct individual hmm_model instances.
 */
void hmm_class::init_groups(const std::vector<metacsst::config::OrderedKeyValues>& motif_groups){
  vector<vector<string>> argv_groups;
  argv_groups.reserve(motif_groups.size());
  for (const auto& group : motif_groups) {
    vector<string> argv;
    argv.push_back("hmm");
    for (const auto& kv : group) {
      const string arg = arg_name(kv.first);
      if (!arg.empty()) {
        argv.push_back(arg);
        argv.push_back(kv.second);
      }
    }
    if (argv.size() > 1) {
      argv_groups.push_back(std::move(argv));
    }
  }

  _number = static_cast<int>(argv_groups.size());
  hmm.clear();
  hmm.reserve(_number);

  for(int i=0;i<_number;i++){
    hmm.push_back(build_hmm(argv_groups[i]));
  }
}

/*build class of HMM models according to the input config file*/
void hmm_class::init(const string& config){
  try {
    const auto dgr_cfg = metacsst::config::parse_dgr_motif_groups(config);
    std::vector<metacsst::config::OrderedKeyValues> merged;
    for (const char* key : {"TR", "VR", "RT"}) {
      const auto it = dgr_cfg.find(key);
      if (it != dgr_cfg.end()) {
        merged.insert(merged.end(), it->second.begin(), it->second.end());
      }
    }
    if (merged.empty()) {
      throw std::runtime_error("No motif groups found in config: " + config);
    }
    init_groups(merged);
    return;
  } catch (const std::exception&) {
  }

  const auto motif_groups = metacsst::config::parse_motif_config(config);
  init_groups(motif_groups);
}

void hmm_class::print(const std::string& dir){
  for(int i=0;i<_number;i++)
    hmm[i].print(dir);

  const auto align_path = std::filesystem::path(dir) / "align.txt";
  const auto score_path = std::filesystem::path(dir) / "score.txt";
  std::ofstream fp1(align_path, std::ios::app);
  std::ofstream fp2(score_path, std::ios::app);
  if(fp1)
    fp1 << "######################################################\n";
  if(fp2)
    fp2 << "######################################################\n";
}

/**
 * @brief Scan a sequence using all clustered GHMMs and merge overlapping hits.
 *
 * Workflow:
 *   1) Run every hmm_model on the target sequence.
 *   2) Collect all match_state hits and sort by start coordinate.
 *   3) Merge overlapping matches on the same strand (scores are summed).
 *   4) Return the unified result set.
 */
std::unique_ptr<out_result> hmm_class::scan_seq(std::string_view seq){

/*WorkFlow:
1>Every GHMM model is used to scan a new sequence,and reserve all the results
2>Filter the results,and if two matchSeqs overlap,merge the two matchSeqs(add the score,it means a stronger information)
3>putout the merged results
*/

  std::vector<match_state> matches;
  matches.reserve(S);

  /*scaning for all the HMM models and reserve all the result_tmps*/
  for(int i=0;i<_number;i++){
    auto result_tmp_sub = hmm[i].scan_seq_full(seq);

    if(result_tmp_sub->number != 0)
      for(int j=0;j<result_tmp_sub->number;j++){
        matches.push_back({
          result_tmp_sub->start[j],
          result_tmp_sub->end[j],
          result_tmp_sub->score[j],
          result_tmp_sub->string[j]
        });
      }
  }

  std::sort(matches.begin(), matches.end(), [](const match_state& lhs, const match_state& rhs) {
    return lhs.start < rhs.start;
  });

  auto result = std::make_unique<out_result>();
  result->number = 0;

  int pos = 0;
  for (const auto& match : matches) {
    if(pos == 0){
      result->start[pos] = match.start;
      result->end[pos] = match.end;
      result->score[pos] = match.score;
      result->string[pos] = match.strand;
      pos++;
      result->number ++;
    }
    else if(match.start < result->end[pos-1] && match.strand==result->string[pos-1]){
      //overlap and merge
      result->end[pos-1] = (result->end[pos-1] > match.end)?result->end[pos-1]:match.end;
      result->score[pos-1] += match.score;
    }
    else{
      result->start[pos] = match.start;
      result->end[pos] = match.end;
      result->score[pos] = match.score;
      result->string[pos] = match.strand;
      pos++;
      result->number ++;
    }
  }

  return result;
}


//struct OUT *searchVR(char *seq,struct OUT **TR,int misMatch);

//=============================================================================================
// Full DGR scanning model (scan_model)
//   - Combines three hmm_class instances (TR, VR, RT) to detect complete DGRs
//=============================================================================================

class scan_model{ //main HMM model used to scan the unknown sequence
 private:
  array<hmm_class,3> state; //three sub state:TR/VR/RT
  int gap; //gap between sub HMMs
 public:
  scan_model(): gap(0) {}
  scan_model(const hmm_class& init_TR,const hmm_class& init_VR,const hmm_class& init_RT,int init_gap){
    init(init_TR,init_VR,init_RT,init_gap);
  }

  void init(hmm_class init_TR,hmm_class init_VR,hmm_class init_RT,int init_gap);
  void print(const std::string& dir);
  std::unique_ptr<out_result> scan_seq(std::string_view seq);

  std::unique_ptr<out_result> scanSeq(std::string_view seq) { return scan_seq(seq); }
};

void scan_model::init(hmm_class init_TR,hmm_class init_VR,hmm_class init_RT,int init_gap){
  state[0]=init_TR;state[1]=init_VR;state[2]=init_RT;
  gap=init_gap;
}

void scan_model::print(const std::string& dir){
  for(int i=0;i<3;i++)
    state[i].print(dir);

  const auto score_path = std::filesystem::path(dir) / "score.txt";
  std::ofstream fp2(score_path, std::ios::app);
  if(fp2){
    fp2 << "Gap Length:" << gap << '\n';
  }
}

/**
 * @brief Detect complete Diversity-Generating Retroelements (DGRs).
 *
 * Workflow:
 *   1) Scan for TR (Template Repeat).
 *   2) If TR hits exist, scan for RT (Reverse Transcriptase).
 *   3) If RT hits exist, scan for VR (Variable Repeat).
 *   4) Collect and sort all sub-hits by genomic coordinate.
 *   5) Group hits into search spaces separated by gaps larger than @c gap.
 *   6) Verify that each search space contains at least TR + VR
 *      (type product divisible by 10, since 2*5=10). Mark index=1 when found.
 */
std::unique_ptr<out_result> scan_model::scan_seq(std::string_view seq){
  auto result = std::make_unique<out_result>();
  result->number = 0;result->index = 0;result->total_score = 0.0f;
  //the final result


/*The scaning order is very important for the efficiency.
  Based on the test result, the scaning order is : TR->RT->VR
*/

  std::array<std::unique_ptr<out_result>, 3> scan_sub;

  scan_sub[0] = state[0].scan_seq(seq);
  if(scan_sub[0]->number >0){ //scaning for TR firstly
    scan_sub[2] = state[2].scan_seq(seq);

    if(scan_sub[2]->number >0){ //scaning for RT secondly

      scan_sub[1] = state[1].scan_seq(seq);
      std::vector<int> start_sub;
      std::vector<int> end_sub;
      std::vector<int> type_sub;
      std::vector<float> score_sub;
      std::vector<int> string_sub;
      start_sub.reserve(S);
      end_sub.reserve(S);
      type_sub.reserve(S);
      score_sub.reserve(S);
      string_sub.reserve(S);

      //fetch all the sub matchSequences to some tmp arrayes
      for(int k=0;k<3;k++)
        for(int i=0;i<scan_sub[k]->number;i++){
          start_sub.push_back(scan_sub[k]->start[i]);
          end_sub.push_back(scan_sub[k]->end[i]);
          score_sub.push_back(scan_sub[k]->score[i]);
          string_sub.push_back(scan_sub[k]->string[i]);
          type_sub.push_back(k+1);
        }

      const int number_sub = static_cast<int>(start_sub.size());

      //sort the array according to the start position
      for(int i=0;i<number_sub-1;i++){
        int pos = i,current_start=start_sub[i];
        for(int j=i+1;j<number_sub;j++)
          //compare and set down the exchange position
          if(start_sub[j] < current_start){
            pos = j;
            current_start=start_sub[j];
          }
        if(pos != i){ //exchange the two matchSeqs
          int tmp1=start_sub[i],tmp2=end_sub[i],tmp3=type_sub[i];
          float tmp4=score_sub[i];int tmp5=string_sub[i];
          start_sub[i]=start_sub[pos];end_sub[i]=end_sub[pos];type_sub[i]=type_sub[pos];score_sub[i]=score_sub[pos];string_sub[i]=string_sub[pos];
          start_sub[pos]=tmp1;end_sub[pos]=tmp2;type_sub[pos]=tmp3;score_sub[pos]=tmp4;string_sub[pos]=tmp5;
        }
      }

      int pos = 0;
      for(int i=0;i<number_sub;){
        if(result->number == 0){
          result->start[pos] = start_sub[i];
          result->end[pos] = end_sub[i];
          result->type[pos] = type_sub[i];
          result->score[pos] = score_sub[i];
          result->string[pos] = string_sub[i];
          result->number ++;
          pos++;
          i++;
        }
        else if(start_sub[i]-result->start[pos-1] <= gap){
          //If a TR and a VR overlaps,merge it to a TR.
          if(start_sub[i] < result->end[pos-1] && (type_sub[i]*result->type[pos-1])%3 != 0){
            result->end[pos-1] = (result->end[pos-1]>end_sub[i])?result->end[pos-1]:end_sub[i];
            result->type[pos-1] = 1;
            result->score[pos-1] += score_sub[i];
            i++;
          }
          else{
            result->start[pos] = start_sub[i];
            result->end[pos] = end_sub[i];
            result->type[pos] = type_sub[i];
            result->score[pos] = score_sub[i];
            result->string[pos] = string_sub[i];
            result->number ++;
            pos++;
            i++;
          }
        }
        else{
          //gap length > gap,check for the last search space,whether a DGR can be found.
          int product = 1; //index for the result,if product%30 == 0 at last,there is a intact DGR structure (30=2*3*5)
          result->total_start = result->start[0];
          result->total_end = result->end[0];
          result->total_score = 0.0f;
          for(int j=0;j<result->number;j++){
            result->total_score += result->score[j];
            if(result->end[j] > result->total_end)
              result->total_end = result->end[j];
            if(result->type[j] == 1)
              product *= 2;
            else if(result->type[j] == 2)
              product *= 3;
            else if(result->type[j] == 3)
              product *= 5;
          }
          if(product%10 == 0){ //a DGR is found:TR+VR,maybe no VR
            result->index = 1;
            break;
          }
          else{
            result->number = 0;
            pos = 0;
          }
        }
      }

      if(result->index == 0){
        //no DGR found for the last search space,check the current search space
        int product = 1;
        result->total_start = result->start[0];
        result->total_end = result->end[0];
        result->total_score = 0.0f;
        for(int j=0;j<result->number;j++){
          result->total_score += result->score[j];
          if(result->end[j] > result->total_end)
            result->total_end = result->end[j];
          if(result->type[j] == 1)
            product *= 2;
          else if(result->type[j] == 2)
            product *= 3;
          else if(result->type[j] == 3)
            product *= 5;
        }
        if(product%10 == 0)
          result->index = 1;
      }
    }
  }
  return result;
}

//=============================================================================================
// Backward-compatible type aliases
//=============================================================================================

using MatchState = match_state;
using OUT = out_result;
using HMM = hmm_model;
using HMM_class = hmm_class;
using SCAN = scan_model;

#endif
