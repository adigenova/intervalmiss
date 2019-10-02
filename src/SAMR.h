//
// Created by Alex Di Genova on 25/07/2016.
//

#ifndef LIGER_SAMR_H
#define LIGER_SAMR_H

//STL classes
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

using namespace std;
//to read an orientation of each SAM file
typedef  struct shortread{
    string name="";
    //string gap="";
    string contig="";
    int pos=-1;
    int len=0;
    bool ori=0;
} shortr;

//structure for building and intervalTREE use 80 bites per record
typedef  struct fragmentPE{
    uint32_t ctgid=0;
    uint32_t pos=0;
    uint16_t d=0;
} fragmentPE;

//faster SAM reader
class SAMR {
private:
  //file descriptor for
  int fd;
  //buffer
  char*  buf;//local buffer for file reading
  //default buffer size
  uint BUFFER_SIZE = 1024*1024*100;//~100Mb real	2m35.414s
  //uint BUFFER_SIZE = 1024;//~20Mb real	2m35.414s
  //char vectors to handle the read from the files
  char liner[2000];//maximal line size
  char remain[2000];//maximal line size
  bool is_empty=false;
  //pointer to move in the buffer
  char *init;//point to the beg
  char *pos;//point to current position in the buffer
  size_t bytes_read;//class variable to check how many bites were reads
  uint32_t consumed=0;
  bool remainb=false;

public:
  //default constructor
  SAMR(string file);
  ~SAMR();
  shortread parse_entry(string read);
  //get the next read
  shortread get_next_read();
  //method that replies if still the file as reads to read
  bool sam_has_reads(){return is_empty;};
  //for error messages from file manipulation
  void handle_error(const char* msg);
  //funtion that load the buffer
  void load_buffer();

};

#endif //LIGER_SAMR_H
