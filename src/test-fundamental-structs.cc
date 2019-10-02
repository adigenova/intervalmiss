
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <thread>
#include <random>
#include <iostream>

using namespace std;

//structure for building and intervalTREE use 80 bites per record
typedef  struct fragmentPE{
    uint32_t ctgid=0;
    uint32_t pos=0;
    uint16_t d=0;
} fragmentPE;

//Low Quality Interval LQI
typedef struct LQI{
    int ctgid=-1;
    int  start=-1;
    int  stop=-1;
    int type=-1;//0 header, 1=interior,2=tail
} LowQualityInterval;



void test_sorted_vector(){
//we create 1 million contigs
std::random_device rd; // obtain a random number from hardware
std::mt19937 eng(rd()); // seed the generator
std::uniform_int_distribution<uint32_t> distr(500, 100000); // define the range
std::uniform_int_distribution<uint32_t> insd(300, 700); // define the of fragments lengths distributions

std::unordered_map<uint32_t,uint32_t> ctglen;
uint32_t tlen=0;
uint32_t number_contigs=1000000;
uint32_t maxcontig=0;
for(uint32_t i=0;i<number_contigs; ++i){
  	//auto l=distr(eng);
    ctglen[i]=distr(eng);
    //std::cout << i<<" "<<ctglen[i]<<std::endl;
    tlen+=ctglen[i];
    //we save the largest contig
    if(maxcontig < ctglen[i]){
      maxcontig=ctglen[i];
    }
}
std::cout << "Total contig length "<< tlen<<std::endl;

//for each contigs we create fragments
//we define the total number of fragmments
uint32_t total_frags=20000000;
//fragments are proportional to the ctg len
vector<fragmentPE> fragments;
auto total_fragments=0;
float acum=0;


for(uint32_t i=0;i<number_contigs; ++i){
  auto p=((double)ctglen[i]/(double)tlen) * total_frags;
  //number of fragments for this contig
  auto numfrag=int(p)+1;
  total_fragments+=numfrag;
  //cout << ctglen[i]<<" "<<i <<" "<<p<<" "<<int(p)+1<<endl;
  std::uniform_int_distribution<uint32_t> distr(0, ctglen[i]-400); // define the range
    for(auto j=0; j< numfrag; j++){
        fragmentPE f;
        f.ctgid=i;
        f.pos=distr(eng);
        f.d=insd(eng);
        fragments.push_back(f);
      }
  }
  cout << "Total number of fragments to create "<< total_fragments<<endl;
  cout << "number of fragmetns stored "<<fragments.size()<<endl;

  //we sort the structure by contig and ctg position within the contigs
  sort( fragments.begin( ), fragments.end( ), [ ]( const fragmentPE a, const fragmentPE b)
  {
     //return a.ctgid < b.ctgid;
     return tie(a.ctgid,a.pos) < tie(b.ctgid,b.pos);
  });
  //we populate the contig to pos hash
  unordered_map<uint32_t,uint32_t> ctg2pos;
  for(uint32_t k=0; k<fragments.size(); ++k){
        //cout <<k<<" "<< fragments[k].ctgid<<" "<<fragments[k].pos<<" "<<fragments[k].d<<" "<<ctglen[fragments[k].ctgid]<<endl;
        if(ctg2pos.count(fragments[k].ctgid)==0){
          ctg2pos[fragments[k].ctgid]=k;
        }
  }


/*
for(uint32_t k=0; k<1000; ++k){
    cout << k<< " " <<ctg2pos[k]<<" "<<ctglen[k]<<" "<<fragments[ctg2pos[k]].d <<endl;
}*/

//I have to determine the LQI intervals of bad quality
//we check the time spent in db creation
clock_t begin = clock();
//we init the static contig to maxcontig
vector<int> cov(maxcontig,0);
uint32_t counter=0;
vector<LQI> results;
for(auto x:ctg2pos){
  //we init the vector of positions
     //vector<int> cov(ctglen[x.first],0);
    for(uint32_t p=x.second; p<fragments.size(); p++){
       //auto tmpf=&fragments[p]; //pointer to takes 76.6377 secs
       //auto tmpf=fragments[p];//copy of fragments[p] takes 71.7188 secs
       //just fragments[p] takes 68.3123 secs
       //there is no big difference between the ways that I got a reference or value from the container
        if(fragments[p].ctgid == x.first){
          //we fill the array
          int max=fragments[p].pos+fragments[p].d;
          if(max>=ctglen[x.first] ){
            max=ctglen[x.first]-1;
          }
          //we fill the fragment
          /*for(auto fi=fragments[p].pos; fi<=max; ++fi){ //takes 131.343
              cov[fi]++;
          }*/
            //we increment the fragment pos
          for (auto it=cov.begin()+fragments[p].pos, end=cov.begin()+max; it!=end; ++it) {*it+=1;}//takes 71.5984 secs

        }else{
          break;
        }
    }


    //variables to determine LQI segments
    uint32_t start=0;
    uint32_t stop=0;
    uint32_t current=0;
    uint32_t maxpos=ctglen[x.first];
    //I have to check the invervals
    for(uint32_t i=0;i<maxpos; ++i){
          //cout << x.first <<" "<<i<<" "<<cov[i]<<endl;
          if(cov[i] < 2){

            if(current==0){
              current = i;
              start = i;
            }
            //allow a max dist of 10
            if(i < current+10){
              current=i;
              stop=i;
            }else{
              //is a new interval
              //cout << "LQI " << current_contig<<" "<<start <<" "<<stop<<" "<<stop-start<<endl;
              LQI r;
              r.ctgid=x.first;
              r.start=start;
              r.stop=stop;
              //the mis assembly is at the interior of the contig length
              r.type=2;//means internal type
              //header type
              if(start <= 1000 and stop <= 1000){
                  r.type=0;
              }
              //tail type
              if(stop > maxpos-1000 and start > maxpos-1000){
                  r.type=1;
              }
              //we store the current LQI in the array
              results.push_back(r);
              current=0;
              start=i;
              stop=i;
            }
        }
    }
    //we test if there is the  last interval is  LQI
    if(current){
              //cout << "LQI " << current_contig<<" "<<start <<" "<<stop<<" "<<stop-start<<endl;
              LQI r;
              r.ctgid=x.first;
              r.start=start;
              r.stop=stop;
              //the mis assembly is at the interior of the contig length
              r.type=2;//means internal type
              //header type
              if(start <= 1000 and stop <= 1000){
                  r.type=0;
              }
              //tail type
              if(stop > maxpos-1000 and start > maxpos-1000){
                  r.type=1;
              }
              //we store the current LQI in the array
              results.push_back(r);
    }




    //we print the results just for checking at this stage
    /*for(auto r:results){
      cout << "LQI"<<r.ctgid<<" "<<r.start<<" "<<r.stop<<" "<<r.type<<endl;
    }*/

    //I clean up the static array
    //std::fill(cov.begin(), cov.end(), 0);
    std::fill(cov.begin(),cov.begin()+ctglen[x.first],0);//takes 95 seconds
    //using memset function
    //memset(&cov[0], 0, ctglen[x.first] * sizeof cov[0]); //takes 93.5287 secs
    counter++;
    if(counter%50000 == 0){
      cout << counter << endl;
    }
    //cov.erase(cov.begin(),cov.end());
}

clock_t end = clock();
double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cout << "Time spent in interval check " << elapsed_secs <<" secs "<<endl;
cout << "Total LQI intervals detected "<< results.size()<<endl;

/*
for(auto r:results){
  cout << "LQI"<<r.ctgid<<" "<<r.start<<" "<<r.stop<<" "<<r.type<<endl;
}*/

}

int main (int argc, char* argv[]){
  test_sorted_vector();
}
