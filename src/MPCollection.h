//
// Created by Alex Digenova on 2/21/18.
//

#ifndef HSCAFF_DEV_MPCOLLECTION_H
#define HSCAFF_DEV_MPCOLLECTION_H


#include <vector>
#include <algorithm>
#include <thread>


class MPLibrary;
#include "MPLibrary.h"
#include "Contig.h"
#include "Finterval.h"
#include "SAMR.h"


using namespace std;

class MPCollection {
private:
    vector<MPLibrary*> libs;
    string filelibs;
    bool is_sort;

public:
    MPCollection(string file);
    void sort_libs();
    void read_libs(Contig *a,Finterval *treetmp,int ncpu);
    void load_lib(int lid, Contig *a);
    void print_frags_libs();
    //void compute_ctg_coverage(Contig *a);
    vector<MPLibrary *> get_all_libs();
    int get_std_by_lib(int libid){return this->libs[libid]->getStd_insert_size();};
    int get_avg_by_lib(int libid){return this->libs[libid]->getAvg_insert_size();};
    void scale_estimated_std(float factor);
};

#endif //HSCAFF_DEV_MPCOLLECTION_H
