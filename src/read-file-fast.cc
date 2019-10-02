#include <algorithm>
#include <iostream>
#include <cstring>

// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "SAMR.h"
//we change to the std
using namespace std;

const char* map_file(const char* fname, size_t& length);
const uintmax_t wc(char const *fname);

int main(int argc, char* argv[])
{

  //std::cout <<  wc(argv[1])<<endl;

  auto samreader = new SAMR(argv[1]);
  uint number_reads=0;
  //cout << samreader->sam_has_reads()<<endl;
  shortread f;
  shortread r;
  while(!samreader->sam_has_reads()){
       auto fwd=samreader->get_next_read();
      if(fwd.pos!=-1){
        number_reads++;
        auto rev=samreader->get_next_read();
        number_reads++;
      /*cout  << fwd.contig<<" "<<rev.contig<<" "<<fwd.pos<<" "<<rev.pos
              <<" "<<fwd.ori<<" "<<rev.ori<<" "<<fwd.len<<" "<<rev.len
              <<" "<<fwd.name<<" "<<rev.name<<endl;*/
              //cout<<" "<<fwd.name<<" "<<rev.name<<endl;
              f=fwd;
              r=rev;
      }
      //cout <<number_reads<<" "<<samreader->sam_has_reads()<<endl;
  }

  cout <<"Number of reads in SAMCLASS "<<number_reads<<endl;
  cout  << f.contig<<" "<<r.contig<<" "<<f.pos<<" "<<r.pos
          <<" "<<f.ori<<" "<<r.ori<<" "<<f.len<<" "<<r.len
          <<" "<<f.name<<" "<<r.name<<endl;

}

void handle_error(const char* msg) {
    perror(msg);
    exit(255);
}

const char* map_file(const char* fname, size_t& length)
{
int fd = open(fname, O_RDONLY);
   if (fd == -1)
       handle_error("open");

   // obtain file size
   struct stat sb;
   if (fstat(fd, &sb) == -1)
       handle_error("fstat");

   length = sb.st_size;

   const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
   if (addr == MAP_FAILED)
       handle_error("mmap");

   // TODO close fd at some point in time, call munmap(...)
   return addr;
}


const uintmax_t wc(char const *fname)
{
     // For a modern 7200 RPM drive it ranges from about 80-160 MB/s today.
    //1Mb buffer it could be just 500Mb
    //there is a small diferrence between using a big or a short buffer
    //auto BUFFER_SIZE = 1024*1024*100;//~100Mb 2m27.061
    //auto BUFFER_SIZE = 1024*1024*10;//~10Mb real	2m35.414s
    auto BUFFER_SIZE = 1024*1024*20;//~20Mb real	2m35.414s

    int fd = open(fname, O_RDONLY);
    if(fd == -1)
        handle_error("open");

        // obtain file size
        struct stat sb;
        if (fstat(fd, &sb) == -1)
            handle_error("fstat");

        //length = sb.st_size;
      if(sb.st_size < BUFFER_SIZE){
        //std::cout << sb.st_size<<endl;
        BUFFER_SIZE=sb.st_size-1;
        std::cout << sb.st_size<<" "<<BUFFER_SIZE<<std::endl;
      }
    /* Advise the kernel of our access pattern.  */
    //posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
    //we allocate dinamically the buffer size according to the file size
    char*  buf = new char[BUFFER_SIZE+1];
    uintmax_t lines = 0;

    string filen(fname);
    auto parser=new SAMR(filen);

    char liner[1000];//maximal line size
    char remain[1000];//maximal line size
    uint32_t consumed=0;
    bool remainb=false;

    while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
    {
        if(bytes_read == (size_t)-1)
            handle_error("read failed");
        if (!bytes_read)
            break;
        //useful function to locate character in the buffer
        //menchr
        char* init=buf;//first position of the buffer

        //p is the advanced position
        char *p;
        consumed=0;
        //string tmp;
        for( p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p){
          memcpy(liner,init,p - init);
          liner[p - init]='\0';
          consumed+=p-init;
          std::string tmp(liner);

          if(remainb){
          string tmp2(remain);
          tmp=tmp2+tmp;
          //consumed=0;
          remainb=false;
          }

          init=p+1;

          //cout << tmp.length()<<" "<<tmp[0]<<endl;
          if(tmp.find("@")){
          auto fwd=parser->parse_entry(tmp);
          ++lines;
         }

            tmp.clear();
            //liner[0]='\0';
          }

          if(bytes_read - consumed > 0){
            //there is a remainig in the buffer
            memcpy(remain,init,buf+bytes_read - init);
            remain[buf+bytes_read - init]='\0';
            //std::cout <<"Remain again "<<bytes_read<<" "<<consumed<<" "<<remain<<endl;
            remainb=true;
          }

    }
    delete[] buf;
    return lines;
}
