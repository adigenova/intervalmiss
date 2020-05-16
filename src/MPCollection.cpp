//
// Created by Alex Digenova on 2/21/18.
//

#include "MPCollection.h"

//threads variables
struct parallellib {
    pthread_mutex_t mut;
    vector<int> *task;
    MPCollection *libs;
    Contig *a;
    //Finterval *treetmp;
};

void * threaded_lib(void* args) {
    //casting of variables
    auto *tmp = (parallellib *) args;
    pthread_mutex_t *mutex = &tmp->mut;
    auto libs = (MPCollection *) tmp->libs;
    auto  contigs=(Contig *) tmp->a;
    //auto  treetmp=(Finterval *) tmp->treetmp;
    auto tasks = (vector<int> *) tmp->task;
    bool dojob = 1;//
    int lid = -1;
    while (dojob) {
        //we try to get a file to process
        pthread_mutex_lock(mutex);
        if (tasks->size() > 0) {
            lid = tasks->back();
            tasks->pop_back();
        } else {
            dojob = 0;
        }
        pthread_mutex_unlock(mutex);
        if (dojob) {
            libs->load_lib(lid,contigs);
        }
    }
    return NULL;
}



//create a set of MP libraries and compute the average insert size if necessary
MPCollection::MPCollection(string files) {
    ifstream readFile(files);
    string file = "";
    filelibs=files;
    MPLibrary *lib;
    //foreach lib file we create a MP library
    while (getline(readFile, file)) {
            //we ask if the file has 3 columns
            istringstream tab(file);
            vector<string> results;
            string val;
            //we parse the lib file
            while(getline(tab,val,'\t')) {
                results.push_back(val);
            }

            switch(results.size()){
                case 2:
                    //lib = new MPLibrary(results[0],atoi(results[1].c_str()),int(atoi(results[1].c_str())*0.1));
                    //wait a type 0 short != long
                    lib = new MPLibrary(results[0],atoi(results[1].c_str()));
                    break;
                case 3:
                    lib = new MPLibrary(results[0],atoi(results[1].c_str()),atoi(results[2].c_str()));
                    break;
                default:
                    lib=new MPLibrary(file);
                    break;
            }
        //we attach the current lib to the Collection
            libs.push_back(lib);
        }

    cout << "Sorting libs by insert size:"<<endl;
    this->sort_libs();
    for(auto l:libs){
        cout << l->getFile() <<" "<<l->getRank() <<" "<<l->getLibType()<<" "<<l->getAvg_insert_size()<<" "<<l->getStd_insert_size()<<endl;
    }
}



void MPCollection::sort_libs() {
    sort( libs.begin( ), libs.end( ), [ ]( const MPLibrary* lhs, const MPLibrary* rhs)
    {
        return lhs->getAvg_insert_size() < rhs->getAvg_insert_size();
    });
    is_sort=true;
    //after sorting we set the rank for each library
    // the order is from smaller to largert insert size
    int i=0;
    for (auto l:libs) {
        l->setRank(i);
        i++;
    }
}



//number of thread used to read the libraries
void MPCollection::read_libs(Contig *a,Finterval *treetmp, int ncpu) {

    if(ncpu == 1) {
      cout <<"reading libraries"<<endl;
        for (auto l:libs) {
            l->load_library(a);
        }
    }else{
        //we load the libs using more than one cooord
        vector<int> cnstask;
        for(int i =0; i < this->libs.size(); i++){
            cnstask.push_back(i);
        }

        pthread_mutex_t mut;
        pthread_mutex_init(&mut, NULL);

        parallellib tmp;
        tmp.mut=mut;//MUTEX to control access to files
        tmp.task=&cnstask;
        tmp.libs=this;
        tmp.a=a;
        //tmp.treetmp=treetmp;

        pthread_t *tab_threads= new pthread_t [ncpu];
        //create the threads
        for(int ii=0;ii<ncpu;ii++) {
            pthread_create(&tab_threads[ii], NULL, threaded_lib, &tmp);
        }
        //wait for thread to finish
        for(int ii=0;ii<ncpu;ii++)
        {
            pthread_join(tab_threads[ii], NULL);
        }
    }

    //we load the treetmp structure and we empty the other positions
    for(auto l:libs){
        auto f=l->get_frags();
        cout << "Lib: "<<l->getRank()<<endl;
        cout << "   Number of pair-ends : "<<l->getNumber_pairs()<<endl;
        cout << "   Number of fragments : "<<f.size()<<endl;

        treetmp->add_fragments(f);
        //we clean the library fragments
        l->clean_frags();
    }


}

void MPCollection::load_lib(int lid, Contig* a) {
    this->libs[lid]->load_library(a);
}

void MPCollection::print_frags_libs() {
    for (auto l:libs) {
        l->print_fragments();
    }
}



vector<MPLibrary *> MPCollection::get_all_libs() {
    return libs;
}

//function to compute std by avg*factor, factor = 0.1
void MPCollection::scale_estimated_std(float factor) {
    for (auto l:libs){
        auto avg=l->getAvg_insert_size();
        l->setStd_insert_size(int(avg*factor+0.5));
        cout << "Scaling std: "<<l->getRank()<<" "<<l->getAvg_insert_size()<<" "<<l->getStd_insert_size()<<endl;
    }
}


