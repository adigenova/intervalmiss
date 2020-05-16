//
// Created by Alex Di Genova on 25/07/2016.
//

#include "SAMR.h"


SAMR::~SAMR(){
//we delete the buffer and close the opened file
	delete[] buf;
	close(fd);

}
//default constructor
SAMR::SAMR(string file){

	fd = open(file.c_str(), O_RDONLY);
	if(fd == -1)
			handle_error("open");
	// obtain file size
	struct stat sb;
	if (fstat(fd, &sb) == -1)
			handle_error("fstat");
	//cout << BUFFER_SIZE<<endl;
	//length = sb.st_size;
	if(sb.st_size < BUFFER_SIZE){
			//std::cout << sb.st_size<<endl;
			BUFFER_SIZE=sb.st_size-1;
			//cout << "ADjusting the buffer size to : "<< sb.st_size<<" "<<BUFFER_SIZE<<endl;
	}
	//we create the buffer
	buf = new char[BUFFER_SIZE+1];
	//we read an initial portion of the file to fill the buf
	bytes_read = read(fd, buf, BUFFER_SIZE);
	//we dont reads lines from the file
	if(bytes_read == (size_t)-1)
			handle_error("read failed");
	//we check that the file is not empty
	if (!bytes_read){
			//break;
			//there is no read in the file
			is_empty=true;
	}
	//we adjust the pointers to the start of the buffer
	init=buf;
	pos=buf;
	consumed=0;

}

//read data from the file
void SAMR::load_buffer(){
  bytes_read = read(fd, buf, BUFFER_SIZE);
	//we dont reads lines from the file
	if(bytes_read == (size_t)-1)
			handle_error("read failed");
	//we check that the file is not empty
	if (!bytes_read){
			//there is no read in the file
			is_empty=true;
			//means that we read data from the file
	}else{
			is_empty=false;
	}
	//we adjust the pointers to the start of the buffer
	init=buf;
	pos=buf;
	//and the consumed bases also
	consumed=0;

}

shortread SAMR::get_next_read(){
	shortread hit;//we create a hit
	//we check that the buffer
	/*if(is_empty){
		//we load the buffer
		load_buffer();
	}*/
	//cout << "in get_next"<<endl;
	//we look for the newline
	pos = (char*) memchr(pos, '\n', (buf + bytes_read) - pos);

	//cout << "in get_next"<<endl;
	//cout << pos<<endl;
	//we found the "\n" in the characters
	if(pos != NULL){
	//we copy the new line to the buffer
	memcpy(liner,init,pos - init);
	//we add the termination character
	liner[pos - init]='\0';
	//we count the consumed based
	consumed+=pos-init;
	//we create the string seq
	string tmp(liner);
  init=pos+1;
	//we add the part of the remain sequence
	/*if(remainb){
		string tmp2(remain);
		tmp=tmp2+tmp;
		remainb=false;
	}*/
	pos++;
	//cout << tmp<<endl;
	if(tmp.find("@"))
	   hit=parse_entry(tmp);

}else{
	//we save the remainder of the buffer and we load the buffer again
	if(bytes_read - consumed > 0){
		//there is a remainig in the buffer
		memcpy(remain,init,buf+bytes_read - init);
		remain[buf+bytes_read - init]='\0';
		//cout <<"Remain again "<<bytes_read<<" "<<consumed<<" "<<remain<<endl;
		remainb=true;
	}
	//we load the buffer for the next read
	load_buffer();
	if(remainb && !sam_has_reads()){
		//we complete the current hit
		pos = (char*) memchr(pos, '\n', (buf + bytes_read) - pos);
		memcpy(liner,init,pos - init);
		//we add the termination character
		liner[pos - init]='\0';
		//we count the consumed based
		consumed+=pos-init;
		//we create the string seq
		string tmp(liner);
	  //init=pos+1;
		init=pos+1;
		//we add the part of the remain sequence
		//if(remainb){
		string tmp2(remain);
		tmp=tmp2+tmp;
		remainb=false;
		pos++;
		//cout << tmp<<endl;
		if(tmp.find("@"))
		   hit=parse_entry(tmp);
	}

}
//we return the seq anyway
	return hit;

}

//printing errors and other things
void SAMR::handle_error(const char* msg) {
    perror(msg);
    exit(255);
}







shortread SAMR::parse_entry(string original){
	shortread rp;
    /*istringstream tab(read);
    vector<string> results;
    string val;
    //we split the SAM by tab
		//it split at 115354560/5m6.952 =~375,747 per sec
    while(getline(tab,val,'\t')) {
        results.push_back(val);
    }*/
		//
		//it split at 115354560/2m56.626 =~651,720 per sec 1.7x faster
		//if we break the while it took 115354560/2m27.061s=~2X faster
		char separator='\t';
		std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    std::string::const_iterator next = std::find( start, end, separator);
		uint32_t counter=0;
    while ( next != end ) {
        results.push_back( std::string( start, next ) );
        start = next + 1;
        next = std::find( start, end, separator );
				counter++;
				//currently we are using the first 10 colum of the SAM
				if(counter > 10){
					break;
				}
    }
    //results.push_back( std::string( start, next ) );
		/*vector<string> results;
		string delimiters = "\t";
		size_t current;
		size_t next = -1;
		do
		{
		  current = next + 1;
		  next = original.find_first_of( delimiters, current );
		  //cout << s.substr( current, next - current ) << endl;
			results.push_back( original.substr( current, next - current ) );
		}
		while (next != string::npos);*/
		//it split at 115354560/3m14.943s =~594,611 per sec 1.7x faster
		/*char separator='\t';
		auto chars_array = strtok(original.c_str(), separator);
		vector<string> results;
    while(chars_array)
    {
				string tmp(chars_array);
				results.push_back(tmp);
        chars_array = strtok(NULL, separator);
    }*/

    //I need an assert to check the the SAM was well parsedl
    assert(results.size()>=6);
    //read_id,orientation,ctg_id,pos,mapq,ksize,qcov
    //parsing each entry in the SAM
    rp.name=results[0];
    rp.contig=results[2];
    //read orientation for single and paired reads aligments
    rp.ori= (atoi(results[1].c_str()) & 16) ? false : true;
    //rp.len=int(results[9].length());
    rp.len=atoi(results[5].c_str());
    rp.pos=atoi(results[3].c_str());
    results.erase(results.begin(),results.end());
	return rp;
}
