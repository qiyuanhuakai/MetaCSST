#ifndef GHMM_MODERN_H
#define GHMM_MODERN_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <iomanip>

///////////////////////////////////////////////////////////////////////////////////////////////
//Package: Metagenomic Complex Sequence Scanning Tool (MetaCSST)                             //
//Developer: Fazhe Yan                                                                       //
//Email: fazheyan33@163.com / ccwei@sjtu.edu                                                 //
//Department: Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University  //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "fun_modern.h"
using namespace std;

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

struct box { //every box is a state in the Hidden Markov Model
  int length; //state matrix length
  int pos_start; //start position of the pattern box in the training set sequence
  int pos_end; //end position of the pattern box in the training set sequence
  float max,min; //max score and min score
  vector<vector<int>> align; //align matrix of the state box
  vector<vector<float>> score; //scoring matrix of the state box
};

struct sub_hmm { //sub HMM strctures ,such as TR/VR/RT
  int start; //start site of the sub HMM in the input sequence
  int end; //end site
  float score; //match score of this sub HMM structure
  int index; //index=-1 -> init;  index=0->TR;  index=1->VR  index=2->RT;
};

struct OUT { //scaning result,for sub HMM(TR/VR/RT) or the total DGR
  int number; //matchSeq number
  float score[S]; //score of the matches
  int start[S]; //start site
  int end[S]; //end site
  int string[S]; //string of the match,1::'+' or 2::'-'

  int type[S]; //used only for scaning fot total DGR;1->TR,2->VR,3->RT;
  int total_start; //used only for scaning for DGR
  int total_end; //used only for scaning for DGR
  float total_score; //used only for scaning for DGR
  int index; //used only for DGR;0->no hot,1->DGR in found
};

class HMM{
 private:
  int len; //length cuttof used to ensure a state
  int size; //state number in the Hidden Markov Model
  int window; //max box length,used as a window size
  int gap; //max gap length between two states
  float cuttof; //cuttof value for a subsequence matching a state box(0~1)
  vector<string> name; //state names
  vector<box> state; //every state is a box,represented by a scoring metrix(Position Specific Scoring Metrix)
  vector<vector<float>> trans; //Transition probability Metrix beteen the states
  vector<float> start,end; //start and end probability for the states
 public:
  float score_cuttof; //score cottof for a sequence belong to this HMM model

  HMM(): len(0), size(0), window(0), gap(0), cuttof(0.0f), score_cuttof(0.0f) {}

  HMM(const vector<vector<float>>& transition,
      const vector<pattern>& metrix,
      int number,
      float value,
      int window_size,
      int gap_length,
      int state_length,
      float seq_score_cuttof){
    init(transition,metrix,number,value,window_size,gap_length,state_length,seq_score_cuttof);
  }

  void init(const vector<vector<float>>& transition,
            const vector<pattern>& metrix,
            int number,
            float cuttof_value,
            int window_size,
            int gap_length,
            int state_length,
            float seq_score_cuttof);
  void print(char *dir);
  struct OUT *scanSeqSingle(char *seq); //only scan for positive string
  struct OUT *scanSeqFull(char *seq); //scan for the both two directions
};

void HMM::init(const vector<vector<float>>& transition,
               const vector<pattern>& metrix,
               int number,
               float cuttof_value,
               int window_size,
               int gap_length,
               int state_length,
               float seq_score_cuttof){
  len = state_length;
  gap = gap_length;
  size = number;
  window = window_size;
  cuttof = cuttof_value;
  score_cuttof = seq_score_cuttof;

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

void HMM::print(char *dir){

  char align[1000],score_file[1000];
  sprintf(align,"%s/align.txt",dir);
  sprintf(score_file,"%s/score.txt",dir);

  FILE *fp1 = fopen(align,"a");
  FILE *fp2 = fopen(score_file,"a");
  if(fp1 == NULL || fp2 == NULL){
    if(fp1 != NULL) fclose(fp1);
    if(fp2 != NULL) fclose(fp2);
    return;
  }

  for(int i=0;i<size;i++)
    if(state[i].length >= len){
      fprintf(fp1,"align matrix:\n");
      for(int j=0;j<4;j++){
        switch(j){
        case 0:fprintf(fp1,"A\t");break;
        case 1:fprintf(fp1,"T\t");break;
        case 2:fprintf(fp1,"C\t");break;
        case 3:fprintf(fp1,"G\t");break;
        }

        for(int k=0;k<state[i].length;k++)
          if(k == state[i].length -1)
            fprintf(fp1,"%d\n",state[i].align[j][k]);
          else
            fprintf(fp1,"%d\t",state[i].align[j][k]);
      }
      fprintf(fp2,"scoring matrix:\n");
      for(int j=0;j<4;j++){
        switch(j){
        case 0:fprintf(fp2,"A\t");break;
        case 1:fprintf(fp2,"T\t");break;
        case 2:fprintf(fp2,"C\t");break;
        case 3:fprintf(fp2,"G\t");break;
        }

        for(int k=0;k<state[i].length;k++)
          if(k == state[i].length -1)
            fprintf(fp2,"%0.2f\n",state[i].score[j][k]);
          else
            fprintf(fp2,"%0.2f\t",state[i].score[j][k]);
      }
    }

  fprintf(fp2,"Start Probabilities\t");
  for(int i=0;i<size;i++)
    if(state[i].length >= len)
      fprintf(fp2,"%0.2f\t",start[i]);
  fprintf(fp2,"\nEnding Probabilities\t");
  for(int i=0;i<size;i++)
    if(state[i].length >= len)
      fprintf(fp2,"%0.2f\t",end[i]);
  fprintf(fp2,"\nTrasnition Probability Matrix\n");
  for(int i=0;i<size;i++)
    if(state[i].length >= len){
      for(int j=0;j<size;j++)
        if(state[j].length >= len)
          fprintf(fp2,"%0.2f\t",trans[i][j]);
      fprintf(fp2,"\n");
    }

  fprintf(fp2,"Window size:%d\n",window);
  fprintf(fp2,"Gap length:%d\n",gap);
  fprintf(fp2,"Match score cuttof:%0.2f\n",score_cuttof);
  fprintf(fp1,"++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  fprintf(fp2,"++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  fclose(fp1);fclose(fp2);
}


struct OUT* HMM::scanSeqSingle(char *seq){
  /*workflow of scaning for TR,VR or RT:
    (1)foreach sub state,set down the matching subSeqs(position as well as score).
    (2)according to the gap length,the whole sequence will be splitted to some search space.
    (3)foreach search space,based on the subSeqs,build many paths according to the veterbi algorithm.
    Every path is a solution for the problem,and the path with the best score will be choosed.
    (4)save all the result satisfying the requirements:score(path) > score_cuttof
   */

  struct OUT *result=(struct OUT *)calloc(1,sizeof(struct OUT));
  result->number=0;

  int S2 = static_cast<int>(strlen(seq));

  vector<int> state_pos(S2,-1);
  vector<int> state_index(S2,-1);
  vector<float> state_score(S2,0.0f);
  int number = 0;

  /*Motif scaning strategy:
    if there is a sequence matched,then skip length will be the length of this motif rather then 1.This methos is able to accelerate the process,but may miss the best match in the mean while,which reduces the whole score,leading to a wrong result.
  */

  if(static_cast<int>(strlen(seq)) > window){
    for(int i=0;i<S2-window-1;){ //scan the sequence for every stat,set down the position as well as matching score
      int index=0;
      //index:in this position,exists a state match? 0:no  1:yes
      for(int j=0;j<size;j++){
        if(state[j].length >= len) {
          float tmp=0.0f;

          for(int k=0;k < state[j].length;k++){
            switch(seq[i+k]){
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

          if(tmp > cuttof*state[j].max && tmp > 0){
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
    if(score_search > score_cuttof){
      result->score[num] = score_search;
      result->start[num] = start_search;
      result->end[num] = end_search;
      num++;
      result->number ++;
    }
  }

  return result;
}

struct OUT* HMM::scanSeqFull(char *seq){
  struct OUT *result=(struct OUT *)calloc(1,sizeof(struct OUT));
  result->number=0;

  struct OUT *result1 = scanSeqSingle(seq); //positive chain

  string seq_complementary = metacsst::complementary(string(seq));
  char* seq_comp_cstr = strdup(seq_complementary.c_str());
  struct OUT *result2 = scanSeqSingle(seq_comp_cstr);
  free(seq_comp_cstr);

  int length = static_cast<int>(strlen(seq));
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
  free(result1);free(result2);
  return result;
}

HMM buildHMM(int ARGC,char* ARGV[]){
  int L=0; //sequence length in the multiAlignment result
  float cov=0.9f; //coverage cuttof in every position to make a pattern box
  int box_len_cuttof=7; //pattern box length cuttof
  int max_box_length=0; //max length of the insured patterns
  int pattern_number=0; //pattern box number in total
  float state_score=0.4f; //cuttof of the score(ratio) used to ensure a state
  float ratio=1.0f; //TP value to control when testing
  int gap = 100; //gap between the states
  float ic = 0.5f; //IC value cuttof

  HMM hmm;

  if(ARGC < 2){
    usage(ARGV[0]);
  }
  else{
    int i,j,k;
    string in_path;
    for(i=0;i<ARGC;i++)
      if(strcmp(ARGV[i],"-build") == 0 && i+1<ARGC)
        in_path = ARGV[i+1];
      else if(strcmp(ARGV[i],"-cov") == 0 && i+1<ARGC)
        cov = static_cast<float>(atof(ARGV[i+1]));
      else if (strcmp(ARGV[i],"-len") == 0 && i+1<ARGC)
        box_len_cuttof = atoi(ARGV[i+1]);
      else if (strcmp(ARGV[i],"-score") == 0 && i+1<ARGC)
        state_score = static_cast<float>(atof(ARGV[i+1]));
      else if (strcmp(ARGV[i],"-ratio") == 0 && i+1<ARGC)
        ratio = static_cast<float>(atof(ARGV[i+1]));
      else if (strcmp(ARGV[i],"-gap") == 0 && i+1<ARGC)
        gap = atoi(ARGV[i+1]);
      else if (strcmp(ARGV[i],"-ic") == 0 && i+1<ARGC)
        ic = static_cast<float>(atof(ARGV[i+1]));
      else if(strcmp(ARGV[i],"-h") == 0)
        usage(ARGV[0]);

    if(in_path.empty() || cov <= 0 || cov >1){
      usage(ARGV[0]);
      return hmm;
    }
    else{
      (void)ic;

      ifstream in_first(in_path);
      if(!in_first.is_open()){
        usage(ARGV[0]);
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

      vector<pattern> scan(M);  //store some subPattern in the whole scoring matrix
      for(i=0;i<M;i++){  //initialization of the pattern boxes
        scan[i].length=0;
        scan[i].max=0.0f;
        scan[i].min=0.0f;
        scan[i].pos_start=-1;
        scan[i].pos_end=-1;
        scan[i].score.assign(4,vector<float>(P,0.0f));
        scan[i].matrix.assign(4,vector<int>(P,0));
        scan[i].sum.assign(P,0);
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
            if(scan[i].length >= box_len_cuttof){
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
            if(scan[i].length >= box_len_cuttof){
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
      float score_cuttof=cuttof(&score_train_ptr,line,ratio);

      hmm = HMM(transition,scan,pattern_number,state_score,max_box_length,gap,box_len_cuttof,score_cuttof);
    }
  }
  return hmm;
}

class HMM_class { //clusters of GHMM model
 public:
  vector<HMM> hmm;
  //multi similar GHMM models,which belongs to different classes
  int _number; //number of clusters(models)

  HMM_class(): _number(0) {}
  explicit HMM_class(const string& config): _number(0) { init(config); }
  explicit HMM_class(const char *config): _number(0) { init(string(config)); }

  void init(const string& config); //initialization,based on the config file
  void init(char *config){ init(string(config)); }
  void print(char *dir);
  struct OUT *scanSeq(char *seq); //scaning a new sequence
};

/*build class of HMM models according to the input config file*/
void HMM_class::init(const string& config){
  ifstream CONFIG(config); //config file
  if(!CONFIG.is_open()){
    _number = 0;
    hmm.clear();
    return;
  }

  vector<vector<string>> argv_groups;
  string tmp;

  while(getline(CONFIG,tmp)){
    string tmp_new = chomp(tmp);
    if(!tmp_new.empty() && tmp_new.back() == '\r')
      tmp_new.pop_back();

    if(tmp_new.find("[motif]") != string::npos){
      //[motif] means a new motif,and a new GHMM model will be built
      argv_groups.push_back(vector<string>());
      argv_groups.back().push_back("hmm");
    }
    else if(tmp_new.find('=') != string::npos && !argv_groups.empty()){
      string name = array_split(tmp_new,'=',0);
      string content = array_split(tmp_new,'=',1);
      argv_groups.back().push_back(arg_name(name));
      argv_groups.back().push_back(content);
    }
  }

  _number = static_cast<int>(argv_groups.size());
  hmm.clear();
  hmm.reserve(_number);

  for(int i=0;i<_number;i++){
    vector<char*> argv_ptr(argv_groups[i].size(),nullptr);
    for(size_t j=0;j<argv_groups[i].size();j++)
      argv_ptr[j] = argv_groups[i][j].data();
    hmm.push_back(buildHMM(static_cast<int>(argv_ptr.size()),argv_ptr.data()));
    //foreach set of arguments,build a corresponding HMM model
  }
}

void HMM_class::print(char *dir){
  for(int i=0;i<_number;i++)
    hmm[i].print(dir);

  char align[1000],score_file[1000];
  sprintf(align,"%s/align.txt",dir);
  sprintf(score_file,"%s/score.txt",dir);
  FILE *fp1 = fopen(align,"a");
  FILE *fp2 = fopen(score_file,"a");
  if(fp1 != NULL)
    fprintf(fp1,"######################################################\n");
  if(fp2 != NULL)
    fprintf(fp2,"######################################################\n");
  if(fp1 != NULL) fclose(fp1);
  if(fp2 != NULL) fclose(fp2);
}

/*scan the new sequences using the clusters of GHMMs*/
struct OUT* HMM_class::scanSeq(char *seq){

/*WorkFlow:
1>Every GHMM model is used to scan a new sequence,and reserve all the results
2>Filter the results,and if two matchSeqs overlap,merge the two matchSeqs(add the score,it means a stronger information)
3>putout the merged results
*/

  int arr_start[S],arr_end[S];int arr_string[S];
  float arr_score[S];

  int num=0;
  /*scaning for all the HMM models and reserve all the result_tmps*/
  for(int i=0;i<_number;i++){
    struct OUT *result_tmp_sub = hmm[i].scanSeqFull(seq);

    if(result_tmp_sub->number != 0)
      for(int j=0;j<result_tmp_sub->number;j++){
        arr_start[num] = result_tmp_sub->start[j];
        arr_end[num] = result_tmp_sub->end[j];
        arr_score[num] = result_tmp_sub->score[j];
        arr_string[num] = result_tmp_sub->string[j];
        num ++;
      }
    free(result_tmp_sub);
  }

  /*Sort the array according to the start position,using quick sort*/
  q_sort_state(arr_start,arr_end,arr_score,arr_string,0,num-1);

  struct OUT *result=(struct OUT *)calloc(1,sizeof(struct OUT));
  result->number = 0;

  int pos = 0;
  for(int i=0;i<num;i++){
    if(pos == 0){
      result->start[pos] = arr_start[i];
      result->end[pos] = arr_end[i];
      result->score[pos] = arr_score[i];
      result->string[pos] = arr_string[i];
      pos++;
      result->number ++;
    }
    else if(arr_start[i] < result->end[pos-1] && arr_string[i]==result->string[pos-1]){
      //overlap and merge
      result->end[pos-1] = (result->end[pos-1] > arr_end[i])?result->end[pos-1]:arr_end[i];
      result->score[pos-1] += arr_score[i];
    }
    else{
      result->start[pos] = arr_start[i];
      result->end[pos] = arr_end[i];
      result->score[pos] = arr_score[i];
      result->string[pos] = arr_string[i];
      pos++;
      result->number ++;
    }
  }

  return result;
}


//struct OUT *searchVR(char *seq,struct OUT **TR,int misMatch);

class SCAN{ //main HMM model used to scan the unknown sequence
 private:
  array<HMM_class,3> state; //three sub state:TR/VR/RT
  int gap; //gap between sub HMMs
 public:
  SCAN(): gap(0) {}
  SCAN(const HMM_class& init_TR,const HMM_class& init_VR,const HMM_class& init_RT,int init_gap){
    init(init_TR,init_VR,init_RT,init_gap);
  }

  void init(HMM_class init_TR,HMM_class init_VR,HMM_class init_RT,int init_gap);
  void print(char *dir);
  struct OUT *scanSeq(char *seq);
};

void SCAN::init(HMM_class init_TR,HMM_class init_VR,HMM_class init_RT,int init_gap){
  state[0]=init_TR;state[1]=init_VR;state[2]=init_RT;
  gap=init_gap;
}

void SCAN::print(char *dir){
  for(int i=0;i<3;i++)
    state[i].print(dir);

  char score_file[1000];
  sprintf(score_file,"%s/score.txt",dir);
  FILE *fp2 = fopen(score_file,"a");
  if(fp2 != NULL){
    fprintf(fp2,"Gap Length:%d\n",gap);
    fclose(fp2);
  }
}

/*Workflow for scan the whole DGR:
  1>the sub structures(TR,VR and RT) are found using Motif-GHMM method
  2>split the whole space into some smaller search space according to the distribution of the gap length
  3>recall DGR structure for each search space
*/
struct OUT *SCAN::scanSeq(char *seq){
  struct OUT *result=(struct OUT *)calloc(1,sizeof(struct OUT));
  result->number = 0;result->index = 0;result->total_score = 0.0f;
  //the final result


/*The scaning order is very important for the efficiency.
  Based on the test result, the scaning order is : TR->RT->VR
*/

  struct OUT *scan_sub[3];

  scan_sub[0] = state[0].scanSeq(seq);
  if(scan_sub[0]->number >0){ //scaning for TR firstly
    scan_sub[2] = state[2].scanSeq(seq);

    if(scan_sub[2]->number >0){ //scaning for RT secondly

      scan_sub[1] = state[1].scanSeq(seq);
      int number_sub=0,start_sub[S],end_sub[S],type_sub[S];
      float score_sub[S];int string_sub[S];

      //fetch all the sub matchSequences to some tmp arrayes
      for(int k=0;k<3;k++)
        for(int i=0,j=number_sub;i<scan_sub[k]->number;i++,j++){
          start_sub[j] = scan_sub[k]->start[i];
          end_sub[j] = scan_sub[k]->end[i];
          score_sub[j] = scan_sub[k]->score[i];
          string_sub[j] = scan_sub[k]->string[i];
          type_sub[j] = k+1;
          number_sub += 1;
        }

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
      free(scan_sub[1]);
    }
    free(scan_sub[2]);
  }
  free(scan_sub[0]);
  return result;
}

#endif // GHMM_MODERN_H
