//
// Created by Alex Digenova on 7/4/18.
//

#ifndef INTERVALMIS_FINTERVAL_H
#define INTERVALMIS_FINTERVAL_H

#include "SAMR.h"
#include "Contig.h"
#include "Utils.h"

#include <unordered_map>
#include <algorithm>
#include <vector>
#include <thread>

//the main problem was the unordered_map<uint32_t,vector<fragmentPE>> buffer; container
//because the speed is not the same has a vector.

using namespace std;
//tree structure

//Low Quality Interval LQI
typedef struct LQI{
    int ctgid=-1;
    int  start=-1;
    int  stop=-1;
    int type=-1;//0 header, 1=interior,2=tail
} LowQualityInterval;

class Finterval {
private:
//local in memory structure to process files
vector<fragmentPE> buffer;
//save the location of the contigs within the buffer
unordered_map<uint32_t,uint32_t>  ctg2buffer;
//structure holding the Low Quality Intervals
vector<LQI> results;
//pointer to contigs object
Contig *contigs;
//to compute the average fragment length
uint64_t avg_fragment_length=0;
//
uint32_t avg_ins=0;//comes from the library used to validate the library

int minimum_cov=2;
int minimum_ctg_len=500;
//default minmum_coverage for considere a base valid



public:
    Finterval(Contig *a, int minimum_cov);
    void add_fragment(fragmentPE tmp);
    void add_fragments(vector<fragmentPE> &tmp);
    //dump the fragment database
    void dump_fragments();
    //sort all the fragments using a single core
    void sort_fragments();
    //check all the contigs
    void Check_all_contigs();
    //print the LQI results
    void print_LQI();
    void break_missasm(string prefix);
    void set_ins_from_lib(uint32_t ins){this->avg_ins=ins;};
    int getMinimum_cov() const;
    void setMinimum_cov(int minimum_cov);
    int getMinimum_ctg_len() const;
    void mark_repeats_frags(int clib);
    void setMinimum_ctg_len(int minimum_ctg_len);
};


#endif //INTERVALMIS_FINTERVAL_H
