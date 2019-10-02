//
// Created by Alex Digenova on 7/4/18.
//

#include "Finterval.h"

//default construct
Finterval::Finterval(Contig *a,int minimum_cov){
  //a pointer to the contigs variable
  this->contigs=a;
  this->avg_ins=0;
  this->minimum_cov=minimum_cov;
  this->setMinimum_ctg_len(500);//default is 500bp for contigs
}
//we store a fragmentPE
void Finterval::add_fragment(fragmentPE tmp){
    //the buffer is a construction of the  method

    buffer.push_back(tmp);
}

//maybe is faster than the push_back method
void Finterval::add_fragments(vector<fragmentPE> &tmp){
  //v1.insert( v1.end(), v2.begin(), v2.end() );
  buffer.reserve(buffer.size() + tmp.size());
  buffer.insert(buffer.end(),tmp.begin(),tmp.end());

}


void Finterval::sort_fragments(){
  //we sort each sub-array, the arrays can be filled by multiples threads
      sort( buffer.begin( ), buffer.end( ), [ ]( const fragmentPE a, const fragmentPE b)
       {
          //return a.ctgid < b.ctgid;
          return tie(a.ctgid,a.pos) < tie(b.ctgid,b.pos);
       });
       //we populate the ctg2buffer structure
       for(auto i=0; i<buffer.size(); ++i){
          auto y=buffer[i];
          //we save the leftmost position of the contig in the buffer
         if(ctg2buffer.count(y.ctgid) == 0){
           ctg2buffer[y.ctgid]=i;
         }
       }
}


void Finterval::dump_fragments(){
 uint32_t cc=0;
   for(auto y:buffer){
    cout << cc <<" "<<y.ctgid<<" "<<y.pos<<" "<<y.d<<" "<<endl;
    cc++;
   }

}

void Finterval::Check_all_contigs(){
  //pick the largest contig present in the dataset
  uint32_t maxcontig=0;
  uint32_t ctglen=0;
  for(int i =0; i < this->contigs->get_number_contig(); ++i){
      //cnstask.push_back(i);
      ctglen=this->contigs->getcontig_length(i);
      if(maxcontig < ctglen ){
          maxcontig=ctglen;
      }
  }

  //we create the container that will hold all the queries
  clock_t begin = clock();
  //we init the static contig to maxcontig
  vector<int> cov(maxcontig,0);
  uint32_t counter=0;

  unordered_map<uint32_t,uint32_t> ctg2pos;//save the position in the buffer structure
  uint32_t nfrags=0;
  uint32_t good_bases=0;
  uint32_t total_bases=0;
  //cout <<"Size of buffer"<<buffer.size()<<endl;

  for(int ctgid =0; ctgid < this->contigs->get_number_contig(); ++ctgid){
  //we locate the contig lefmost position in all the array buffers
   if(ctg2buffer.count(ctgid) == 0){
     //we skypt the contig if not present in the ctg2buffer structure
     continue;
   }
  //we populate the fragments
  uint32_t ctgsize=this->contigs->getcontig_length(ctgid);
  total_bases+=ctgsize;
    for(auto i=ctg2buffer[ctgid]; i<buffer.size(); ++i){
          auto f=buffer[i];
          if(f.ctgid == ctgid){
            auto stop=f.pos+f.d;
            //we check that we dont over
            if(stop >= ctgsize ){
              stop=ctgsize-1;
            }
            nfrags++;

            avg_fragment_length+=f.d;

            for (auto it=cov.begin()+f.pos, end=cov.begin()+stop; it!=end; ++it) {*it+=1;}//takes 71.5984 secs
          }else{
              break;
          }
      }
  //variables to determine LQI segments
  uint32_t start=0;
  uint32_t stop=0;
  uint32_t current=0;
  uint32_t maxpos=ctgsize;
  uint32_t limit_type=this->avg_ins+this->avg_ins*0.1*4.5;//boundaries for start/stop misstypes
  //I have to check the invervals
  for(uint32_t i=0;i<maxpos; ++i){
        //cout << x.first <<" "<<i<<" "<<cov[i]<<endl;
        //two or lest than two kmers
        if(cov[i] < minimum_cov){
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
            r.ctgid=ctgid;
            r.start=start;
            r.stop=stop;
            //the mis assembly is at the interior of the contig length
            r.type=2;//means internal type
            //header type
            if(start <= limit_type and stop <= limit_type){
                r.type=0;
            }
            //tail type
            if(stop > maxpos-limit_type and start > maxpos-limit_type){
                r.type=1;
            }
            //entire type
            auto  miscov=(float)(stop-start)/(float)ctgsize;
            //cout <<miscov<<endl;
            if(miscov > 0.8){
                r.type=3;
            }else{
              if(start < 10){
                  r.type=0;
              }
              if(stop > maxpos-10){
                  r.type=1;
                }
            }
            //we store the current LQI in the array
            results.push_back(r);
            current=0;
            start=i;
            stop=i;
          }
      }else{
        //are bases that are supportted by at lest two fragments
          good_bases++;
      }
  }
  //we test if there is the  last interval is  LQI
  if(current){
            //cout << "LQI " << current_contig<<" "<<start <<" "<<stop<<" "<<stop-start<<endl;
            LQI r;
            r.ctgid=ctgid;
            r.start=start;
            r.stop=stop;
            //the mis assembly is at the interior of the contig length
            r.type=2;//means internal type
            //header type
            if(start <= limit_type and stop <= limit_type){
                r.type=0;
            }
            //tail type
            if(stop > maxpos-limit_type and start > maxpos-limit_type){
                r.type=1;
            }

            auto  miscov=(float)(stop-start)/(float)ctgsize;
            //cout <<miscov<<endl;
            if(miscov > 0.8){
                r.type=3;
            }else{
              if(start < 10){
                  r.type=0;
              }
              if(stop > maxpos-10){
                  r.type=1;
                }
            }
            //we store the current LQI in the array
            results.push_back(r);
  }
  //we have to clean the cov array
  std::fill(cov.begin(),cov.begin()+maxpos,0);//takes 95 seconds
  counter++;
  if(counter%50000 == 0){
    cout << counter << endl;
  }

}

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  avg_fragment_length=int(avg_fragment_length/nfrags);
  cout << "Total fragments " << nfrags <<endl;
  cout << "Average fragment length " << avg_fragment_length<<endl;
  cout << "Total Examinated bases " << total_bases<<endl;
  cout << "Total Validated bases by at least "<<minimum_cov<<" fragments " << good_bases<<endl;
  cout << "Percetaje of validated bases  " << ((float)good_bases/(float)total_bases * 100)<<endl;
  cout << "Time spent in interval check " << elapsed_secs <<" secs "<<endl;
  cout << "Total LQI intervals detected "<< results.size()<<endl;

}

void Finterval::print_LQI(){
  for(auto r:results){
    cout << "LQI "<<this->contigs->get_contig_name(r.ctgid)<<" "<<r.ctgid<<" "<<r.start<<" "<<r.stop<<" "<<r.type<<" "<<abs(r.stop-r.start)<<" "<<this->contigs->getcontig_length(r.ctgid)<<endl;
    //cout << "LQI "<<" "<<r.ctgid<<" "<<r.start<<" "<<r.stop<<" "<<r.type<<" "<<abs(r.stop-r.start)<<endl;

  }
}
//produce a fasta and a GFF file describing the miss-assemblies, also mark as repetitive or short the quimeric contigs to avoid them in the PE-scaffoldind step
void Finterval::break_missasm(string prefix){
  cout << "Scanning LQI intervals for missassemblies"<<endl;
  //cout << "Minium contig length  "<<minimum_ctg_len<<" bp to report"<<endl;

    //we save the missassemblies to the following file
  ofstream logmis;
  logmis.open(prefix+".LQIM.txt");
  logmis << "# LQI missassembly types:\n# 0 = LQI at the begining of the contig(first avg_fragment_length+std_fragment_lenght*4)\n# 1 = LQI at the end of the contig (last first avg_fragment_length+std_fragment_lenght*4)\n# 2 = LQI at the interior of the contig\n# 3 = LQI represent a significant part of the contig(> 0.7)\n";
  logmis << "#LQI CTG_NAME CTG_ID LQI(START) LQI(STOP) LQI(TYPE) LQI(LENGTH) CTG_LEN\n";
  unordered_map<uint32_t,uint32_t> quimeric;//save a reference to the first miss assembly event
  //we consider a 10% std then we allow a maximum of avg+std*2.5
  uint32_t max_be_allowed=avg_fragment_length+(avg_fragment_length*0.1)*2.5;//max allowed miss assembly at the beginnig/end
  //cout << max_be_allowed<<endl;
  uint32_t misasm=0;
  uint32_t mstart=0;
  uint32_t mstop=0;
  uint32_t minternal=0;
  uint32_t mwhole=0;
  uint32_t index=0;
  vector<LQI> aerrors;
  for(auto r:results){
    if(r.type == 2 ){
      //we save the lefmost position of the potention
      if(quimeric.count(r.ctgid)==0){
        quimeric[r.ctgid]=index;
      }
      aerrors.push_back(r);
      index++;
      logmis << misasm<<" "<<this->contigs->get_contig_name(r.ctgid)<<" "<<r.ctgid<<" "<<r.start<<" "<<r.stop<<" "<<r.type<<" "<<abs(r.stop-r.start)<<" "<<this->contigs->getcontig_length(r.ctgid)<<endl;
      misasm++;

      minternal++;
    }else{
      if(abs(r.stop-r.start) > max_be_allowed){
        if(quimeric.count(r.ctgid)==0){
          quimeric[r.ctgid]=index;
        }
        aerrors.push_back(r);
        index++;
        logmis << misasm<<" "<<this->contigs->get_contig_name(r.ctgid)<<" "<<r.ctgid<<" "<<r.start<<" "<<r.stop<<" "<<r.type<<" "<<abs(r.stop-r.start)<<" "<<this->contigs->getcontig_length(r.ctgid)<<endl;
        //we mark the type of missassemblies
        if(r.type == 0)
        mstart++;
        if(r.type == 1)
        mstop++;
        if(r.type == 3)
        mwhole++;

        misasm++;
      }
    }
    //index++;
  }

  //we close the file pointer
  cout << "A total of "<<misasm<< " missassemblies were detected from LQI segments"<<endl;
  cout << "MisAssemblies at start = "<< mstart << " at end = "<<mstop<<", internal = "<<minternal<<" and whole contig = "<<mwhole<<endl;

  //we create the fasta of the contigs
  ofstream fasta; //save normal ctgs and quimeric splitted
  fasta.open(prefix+".msplit.fa");
  ofstream qsplit;//save only quimeric contigs useful for post-processing
  qsplit.open(prefix+".qmsplit.fa");

  for(int ctgid =0; ctgid < this->contigs->get_number_contig(); ++ctgid){
        auto ctg=this->contigs->getcontig(ctgid);

      if(quimeric.count(ctgid)==0){
        //we skypt too short contigs
        if(ctg.seq.length() < minimum_ctg_len){
          continue;
        }

          string header=">"+ctg.namectg+" "+to_string(ctg.coverage)+"\n";
          auto seqf=seqformat(60,ctg.seq);
          fasta.write(header.c_str(),header.length());
          fasta.write(seqf.c_str(),seqf.length());
        //fasta << ">"<<ctg.namectg<<" "<<ctg.coverage<<endl;
      }else{
        //we have to cut the sequence
        //and mark it as a repeat contig to avoid using it in scaffodling step
        this->contigs->set_repeat(ctgid);//the contig will be ignored during scaffolding
        this->contigs->set_quimeric(ctgid);//the contig will be removed from the scaffolding ouput
        //we move in the reverse order
        //we save the quimeric intervals of the contig as pairs
        //we create a vector holding the ctglen
        vector<int> cov(ctg.length,0);
        for(auto k=quimeric[ctgid]; k<aerrors.size(); ++k){
              auto m=aerrors[k];
              //cout <<m.ctgid<<" "<<m.start<<" "<<m.stop<<" "<<m.type<<" "<<endl;
              if(m.ctgid !=ctgid){
                break;//we break the for
              }
              //we mask the bad interval of the contig
              for(auto b=m.start; b<=m.stop; ++b){
                 cov[b]++;
              }
        }
        uint32_t gstart=0;
        uint32_t gstop=0;
        uint32_t p=1;
        for(auto b=0; b<cov.size(); ++b){
            //cout << gstart<<" "<<gstop<<" "<<p<<endl;
            if(cov[b]==0){
                //we ignore the base
                gstop++;
            }else{
              //gstop++;
              //we print the good fragment if the seq is larger than 500
              if(gstop >= minimum_ctg_len){
              //cout << ctgid<<" "<<gstart<<" "<<gstop<<" "<<p<<endl;
              string header=">"+ctg.namectg+"_"+to_string(p)+" "+to_string(ctg.coverage)+"\n";
              fasta.write(header.c_str(),header.length());
              qsplit.write(header.c_str(),header.length());
              p++;
              //we print the bad seq
              auto subseq=ctg.seq.substr(gstart,gstop);
              //we format the output seq
              auto seqf=seqformat(60,subseq);
              //we print a nice fasta
              fasta.write(seqf.c_str(),seqf.length());
              qsplit.write(seqf.c_str(),seqf.length());
              }
              gstart=b+1;
              //stop comeback to 0
              gstop=0;
            }
        }

        //we print if gstop is larger than 500
        if(gstop >= minimum_ctg_len){
        //cout << ctgid<<" "<<gstart<<" "<<gstop<<" "<<p<<endl;
        string header=">"+ctg.namectg+"_"+to_string(p)+" "+to_string(ctg.coverage)+"\n";
        fasta.write(header.c_str(),header.length());
        qsplit.write(header.c_str(),header.length());
        p++;

        //we print the bad seq
        auto subseq=ctg.seq.substr(gstart,gstop);
        auto seqf=seqformat(60,subseq);
        fasta.write(seqf.c_str(),seqf.length());
        qsplit.write(seqf.c_str(),seqf.length());
        }

      }
    }
    //we close all the opened files
    logmis.close();
    fasta.close();
    qsplit.close();
    //we clean the results
    results.erase(results.begin(),results.end());

}

int Finterval::getMinimum_cov() const {
    return minimum_cov;
}

void Finterval::setMinimum_cov(int minimum_cov) {
    Finterval::minimum_cov = minimum_cov;
}

int Finterval::getMinimum_ctg_len() const {
    return minimum_ctg_len;
}

void Finterval::setMinimum_ctg_len(int minimum_ctg_len) {
    Finterval::minimum_ctg_len = minimum_ctg_len;
}



void Finterval::mark_repeats_frags(void) {
    //we populate the simple hash
    simplehash32 ctg2frags;
    //we fill the hash with the fragment lenght
    for(auto f:buffer){
        ctg2frags[f.ctgid]+=f.d;
    }
    //we iter the contigs to compute the coverage
    vector<float> cov;
    for (auto con:this->contigs->get_all_contigs()) {
        float c=0;
        if(ctg2frags.count(con.ctgid)){
            c=(float)ctg2frags[con.ctgid]/(float)con.length;
            //cout << "Coverage CTG: "<<con.ctgid <<" "<<c<<" "<<con.length<<" "<<ctg2frags[con.ctgid]<<endl;
        }

        if(c>0) {
            cov.push_back(c);
        }
        //we save the coverage of the contig
        this->contigs->set_ctg_coverage(con.ctgid,c);
    }
    //AVG coverage and STD
    float avg_ctg_cov=compute_average(cov);
    float std_ctg_cov=compute_std(cov,avg_ctg_cov);
    //we set the variables of the contig coverage
    this->contigs->setAvg_ctg_cov(avg_ctg_cov);
    this->contigs->setStd_ctg_cov(std_ctg_cov);
    //todo:improve this when the std_ctg_cov is to large
    //this could be inferrend from the illuumina asssembler
    //float max = avg_ctg_cov + 4*std_ctg_cov;
    auto max = avg_ctg_cov + float(4*std_ctg_cov);
    int repeats=0;
    cout << "Max allowed coverage "<<max<<endl;
    int max_rep_size=0;
    for (auto con:this->contigs->get_all_contigs()) {
        //cout << "NOREPEAT CTG:"<<con.ctgid<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
        if(con.coverage > max){
            this->contigs->set_repeat(con.ctgid);
            cout << "REPEAT CTG:"<<con.ctgid<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
            repeats++;
            if(con.length > max_rep_size){
                max_rep_size=con.length;
            }
        }else{
            //cout << "NONREPEAT CTG:"<<con.ctgid<<" "<<con.length<<" "<<con.repeat<<" "<<con.coverage<<endl;
        }
    }
    cout << "AVG "<<avg_ctg_cov <<" STD "<< std_ctg_cov <<" CTG REP "<<repeats<<" Max rep size "<<max_rep_size<<endl;
}
