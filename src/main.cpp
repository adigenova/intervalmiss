#include <iostream>

//local libraries
#include "Contig.h"
#include "MPCollection.h"
//options
#include "imopts.h"

int main (int argc, char* argv[]){


  	gengetopt_args_info ai;
    if (cmdline_parser (argc, argv, &ai) != 0) {
	//cmdline_parser_print_help ();		
	cout << "Run "<<argv[0]<<" -h to see the list of options. "<<endl;
        exit(1);
    }
 
    string prefix=ai.prefix_arg;//prefix for output
    int cores=ai.cpu_arg;
    //load contigs from fai file
    auto contigs=new Contig(ai.contigs_arg);
    //contigs->print_contigs_file(prefix);
    //creates and reads the libraries
    auto libs=new MPCollection(ai.samlist_arg);
    //object of contigs and min depth coverage for a base
    auto treetmp = new Finterval(contigs,ai.mbc_arg[0]);
    //we ask if the fst was given, currenly  is intend only for minia
    if(ai.fst_given){
        libs->scale_estimated_std(ai.fst_arg);
    }


    //we read the libraries
    libs->read_libs(contigs,treetmp,cores);
    //we check that the links were loaded correctly
    //libs->print_link_libs();
    //we use the contig coverage from short-reads and we use it to mark the repetitive contigs;

    //we check for the parameter of minimum contig length for default if 500bp
    if(ai.mcl_given){
        treetmp->setMinimum_ctg_len(ai.mcl_arg);
    }
    //we set the avg to treetmp
    treetmp->set_ins_from_lib(libs->get_avg_by_lib(0));
    //we sort the fragmments within the three
    treetmp->sort_fragments();
    //treetmp->dump_fragments();
    //we use the populated fragments
    /*for(int i=0; i<contigs->get_number_contig(); i++){
      treetmp->Check_contig(i);
    }
    treetmp->dump_fragments();*/

    //todo: improve the computation of coverage from the fragments
    if(ai.ccoverage_given){
        contigs->mark_repeats(ai.ccoverage_arg);
    }else{
        //we compute the contig coverage and we mark the ones that are repeats from long reads
        //contigs->mark_repeats();
        //we use the fragments to compute the physical coverage of the contigs
        treetmp->mark_repeats_frags();
    }

    for (int i=0; i<ai.mbc_given; ++i) {
        cout << "=========== Checking Contigs ========="<<endl;
        cout <<"Minimum Base Coverage : "<<ai.mbc_arg[i]<<endl;
        cout <<"Minimum Contig length : "<<ai.mcl_arg<<endl;
        treetmp->setMinimum_cov(ai.mbc_arg[i]);//we set the minimum base coverage
    	treetmp->Check_all_contigs();//we have to clean the variables of this function
        string p = prefix + ".MBC" + to_string(ai.mbc_arg[i]);
        treetmp->break_missasm(p);//we write the results, the prefix should change at every step
        cout << "===========  Done    ========="<<endl;
        cout <<endl;
    }
    /*
    //check the contigs in parallel using as many treads
    treetmp->Check_all_contigs();
    //we print the results and quit
    //treetmp->print_LQI();
    //we have to mark quimeric contigs an split them
    treetmp->break_missasm(prefix);
    */

    return 0;
}
