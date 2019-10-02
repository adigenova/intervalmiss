//
// Created by Alex Digenova on 2/20/18.
//

#ifndef HSCAFF_DEV_MPLIBRARY_H
#define HSCAFF_DEV_MPLIBRARY_H

#include <vector>
#include <map>
#include <math.h>
#include<algorithm>


#include "Contig.h"
#include "SAMR.h"
#include "Utils.h"
#include "Finterval.h"

using namespace std;

//struct for links amongs contigs
typedef  struct linking{
    int lid=0; // link id
    //short int lib_id=0;
    /*int ctg1=-1;
    int ctg2=-1;*/

    uint32_t c1=0;
    uint32_t c2=0;
    //link positions
    int p1=-1;
    int p2=-1;
    //inferred contig orientation
    bool cor1=0;
    bool cor2=0;
    //actual mate orientation from alignment
    bool ror1=0;
    bool ror2=0;
    //distance between contigs
    int dist=0;
    bool switched=0;
    linking(int idl, uint32_t c1, uint32_t c2,int p_1, int p_2, bool co1, bool co2, bool ro1, bool ro2, int d,bool sw )
            :lid(idl),c1(c1),c2(c2),p1(p_1),p2(p_2), cor1(co1),cor2(co2),ror1(ro1),ror2(ro2),dist(d),switched(sw) {};
} linkingp;


class MPLibrary {

private:
    int avg_insert_size;
    int std_insert_size;
    int rank;//libraries are ranked according to their average insert size
    int mobs;//store the maximal number of obs to infer average insert size
    int number_pairs;

    float pout; //fractions of outliers to discards while computing the average insert sizes
    string file;//file saving the SAM file
    void _infer_insert_size();
    //save the long reds that pass througth the contig
    //map<string, uint16_t > long_reads;
    vector<linking> links;
    vector<fragmentPE> frags;//container to store locally fragmments

public:
    //default constructor
    MPLibrary(string file);
    //another constructor which specified the outlayers
    MPLibrary(string file, int avg_ins, int std_ins);

    //Getter & Setter methods
    int getAvg_insert_size() const {
        return avg_insert_size;
    }

    void setAvg_insert_size(int avg_insert_size) {
        MPLibrary::avg_insert_size = avg_insert_size;
    }

    int getStd_insert_size() const {
        return std_insert_size;
    }

    void setStd_insert_size(int std_insert_size) {
        MPLibrary::std_insert_size = std_insert_size;
    }

    int getRank() const {
        return rank;
    }

    void setRank(int rank) {
        MPLibrary::rank = rank;
    }

    const string &getFile() const {
        return file;
    }

    void setFile(const string &file) {
        MPLibrary::file = file;
    }

    int getNumber_pairs() const;

    void setNumber_pairs(int number_pairs);

    // void compute_lib_std(vector<int> obs_inserts);
    //function that load the library and discard the pair
    void load_library(Contig *a);

    vector<fragmentPE> get_frags(){return this->frags;};
    void clean_frags(){this->frags.erase(this->frags.begin(),this->frags.end());};
    void print_fragments();
};


#endif //LIGER_MPLIBRARY_H
