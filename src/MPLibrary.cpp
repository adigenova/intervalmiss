//
// Created by Alex Digenova on 2/20/18.
//

#include "MPLibrary.h"


//default construct
MPLibrary::MPLibrary(string file) {
    //we init the defaults variables
    mobs=100000;
    avg_insert_size=0;
    std_insert_size=0;
    rank=0;
    pout=0;
    number_pairs=0;
    this->file=file;
    //we infer the insert size of the library
    this->_infer_insert_size();
}


MPLibrary::MPLibrary(string file, int avg_ins, int std_ins ) {
    mobs=100000;
    avg_insert_size=avg_ins;
    std_insert_size=std_ins;
    rank=0;
    pout=0.1;
    this->file=file;
}


void MPLibrary::_infer_insert_size(){
    auto samreader=new SAMR(this->file);
    vector<int> obs_insert;
    int counter=0;
    while(!samreader->sam_has_reads()){
        auto fwd=samreader->get_next_read();
        //we are reading the header of the SAMfile
        if(fwd.pos==-1) {
            continue;
        }
        //auto fwd=parser.parse_entry(lfwd);
        auto rev=samreader->get_next_read();
        //we are in the same contig
        if(rev.contig.compare(fwd.contig) == 0){
                obs_insert.push_back(abs(rev.pos-fwd.pos));
                counter++;
        }
        //we have all the observations that we need
        if(counter >=mobs){
            break;
        }
    }
    //we close the current SAM file
    //readsam.close();
    //we call the destructor
    samreader->~SAMR();
    //we compute the avg insert size distance
    //cout << "Num obs: "<<counter << endl;
    //this->compute_lib_average(obs_insert);
    this->avg_insert_size=compute_average(obs_insert);
    //this->compute_lib_std(obs_insert);
    //this->std_insert_size=int(avg_insert_size * 0.1);
    this->std_insert_size=compute_std(obs_insert,this->avg_insert_size);
}


void MPLibrary::load_library(Contig *a) {
    auto samreader=new SAMR(this->file);
    //local fragments
    //vector<fragmentPE> localfragments;
    //boundaries for max and min distance among pairs
    float min = float(-this->std_insert_size * 2.5);
    float max = float(this->avg_insert_size+this->std_insert_size * 2.5);
    float minf = float(this->avg_insert_size-this->std_insert_size * 2.5);
    float maxf = float(this->avg_insert_size+this->std_insert_size * 2.5);

    int idl=0;//and id for each stored link

    uint counter=0;
    while(!samreader->sam_has_reads()){
        auto fwd=samreader->get_next_read();
        //we are reading the header of the SAMfile
        if(fwd.pos==-1) {
            continue;
        }
        counter++;
        auto rev=samreader->get_next_read();
        //we are in the same contig
        if(rev.contig.compare(fwd.contig) == 0){
            //we save the contig coverage of the shorter library for repeat identification
            if(this->rank ==0){
                a->add_ctg_reads(a->ctg2id(fwd.contig), 2);
            }
            //both orientations are diferents within contigs valid pairs are FR or RF but never FF or RR
             if(fwd.ori != rev.ori){
                  auto d=abs(fwd.pos-rev.pos);
                  //values for within contigs  reads
                  if(minf < d and maxf > d){
                    fragmentPE ftmp;
                    //we have a valid fragment mapped within the contig
                    //uint16_t pf=0;
                    auto pf=fwd.pos;
                    if(fwd.pos > rev.pos){
                      pf=rev.pos;
                    }
                    //this is the fragment information for building a contig file
                    //cout << fwd.contig <<" "<<pf<<" "<<d<<endl;
                    ftmp.ctgid=a->ctg2id(fwd.contig);
                    ftmp.pos=pf;
                    ftmp.d=uint16_t(d);
                    frags.push_back(ftmp);
                  }
              }
        }
    }

    this->setNumber_pairs(counter);//store the number of pairs read in the library
    //we remove the samreader
    samreader->~SAMR();

}

//funtion to print the links of a given library
void MPLibrary::print_fragments(){

    for(auto l:this->frags){
            cout <<"OK "<<this->getRank()<<" "<<l.ctgid << " " << l.pos <<" "<<l.d<<  endl;
    }

}


int MPLibrary::getNumber_pairs() const {
    return number_pairs;
}

void MPLibrary::setNumber_pairs(int number_pairs) {
    MPLibrary::number_pairs = number_pairs;
}
