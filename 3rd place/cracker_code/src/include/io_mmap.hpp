/*
* using mmap for intensive io
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <vector>

using namespace std;

/* 
 * record the file block 
 * each has 5G.
 * */
class fileblock{
    public:
        long offset;
        long length;
        long shift;

        fileblock():offset(0), length(0), shift(0){}
        fileblock(long offset_, long len_, long shift_):
            offset(offset_), length(len_), shift(shift_){}
};

class mmap_reader{
    public:
        int fd;
        char *start;
        struct stat sb;
        long end; /* offset */
        long length;

        /* record the split file*/
        bool is_split;
        vector<fileblock> fb;
        
        mmap_reader():end(0), is_split(false){}

    /* normally open a file, when it is less than 5G. */
bool open_mmap(char* filepath){
    end = 0;
  fd = open(filepath, O_RDONLY);
  fstat(fd, &sb);
  length = sb.st_size;
//  printf("filesize=%ld\n", sb.st_size);
  start = (char*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
  if(start == MAP_FAILED)
      return false;
  return true;
}

/* open a file block */
bool open_mmap(char* filepath, int blk_idx){
  fd = open(filepath, O_RDONLY);
  fstat(fd, &sb);
  length = fb[blk_idx].length;
  long offset = fb[blk_idx].offset - fb[blk_idx].shift;
  end = fb[blk_idx].shift;
//  printf("filesize=%ld\n", sb.st_size);
  start = (char*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, offset);
  if(start == MAP_FAILED){
     //printf("error number=%d\n", errno);
      return false;
  }
  return true;
}

/* remap a file if it is no split for a file larger than 5G. */
bool remap(char* filepath, long offset){
  munmap(start, length);
  size_t page_size = sysconf(_SC_PAGESIZE);
  long move_back = offset % page_size;
  length = sb.st_size - offset + move_back;

//  printf("page_size=%d\n", page_size);
//  printf("files=%ld, length=%ld, offset=%ld move_back=%ld noffset=%ld\n", sb.st_size, length, offset, move_back, offset-move_back);
  start = (char*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, offset - move_back);
  end = move_back;
  if(start == MAP_FAILED){
 //     printf("error nonumbr=%d\n", errno);
      return false;
  }
  return true;
}

int get_fileblock_num(){
    return fb.size();
}

bool split_file(char* dir){
    long limit = 5L * 1024 * 1024 * 1024;
    fd = open(dir, O_RDONLY);
    fstat(fd, &sb);
    if(sb.st_size <= limit){
        is_split = false;
        close(fd);
        return false;
    }
    is_split = true;
    fb.clear();
    fstat(fd, &sb);
    long length = 0, shift = 0, noffset = 0;
    size_t page_size = sysconf(_SC_PAGESIZE);
    //printf("lim=%ld fs=%ld\n", limit, sb.st_size);
    while(length + noffset < sb.st_size){
       noffset = noffset + length;
       length = limit;
       if(noffset + length > sb.st_size){
          length = sb.st_size - noffset; 
       }
       else{
//           long curoffset = tell(fd);
           //lseek64(fd, noffset + length - curoffset, SEEK_CUR);
           lseek64(fd, noffset + length, SEEK_SET);
           char a;
           int cnt = 0;
           while(read(fd, &a, 1) != 0){
               ++cnt;
               if(a == '\n') break;
           }
           length += cnt;
           //printf("cnt=%d\n", cnt);
       }
       shift = noffset % page_size;
//       printf("offset=%ld len=%ld shift=%ld next=%ld\n", noffset, length, shift, noffset+length);
       fb.push_back(fileblock(noffset, length + shift, shift));
    }
    //printf("fileblock num = %d\n", fb.size());
    close(fd);
    return true;
}

void close_mmap(){
  munmap(start, length);
  close(fd);
}

long filesize(){
    return sb.st_size;
}

long map_file_size(){
    return length;
}

long linecount(){
    long cnt = 0;
    while (end < sb.st_size){
        while(start[end++] !='\n');
        cnt++;
    }
    return cnt;
}
long skipfirstline(){
  end = 0; 
  while (start[end++]!='\n');
  return end;
}

long skipline(){
    long old = end;
  while (start[end++]!='\n');
  return end - old;
}

long skiplinefrom(long& offset){
  long old = offset;
  while (start[offset++]!='\n');
  return offset - old;
}

long skiprecord(){
  long old = end;
  while(start[end] != '|' && start[end] != '\n') 
    ++end;
  ++end;
  return end - old;
}


bool hasNext(){
    return end < length;
}

/* read a integer */
int getInt() {
  int x = 0;
  while (start[end]>='0' && start[end]<='9') {
    x = x*10+(start[end]-'0');
    ++end;
  }
  ++end; /* skip '|' */
  return x;
}

long getInt(int& x) {
  long old = end;
  x = 0;
  while (start[end]>='0' && start[end]<='9') {
    x = x*10+(start[end]-'0');
    ++end;
  }
  ++end; /* skip '|' */
  return end - old;
}

long getLong() {
  long x = 0;
  while (start[end]>='0' && start[end]<='9') {
    x = x*10+(start[end]-'0');
    ++end;
  }
  ++end; /* skip '|' */
  return x;
}

long getLong(long& x) {
  long old = end;
  x = 0;
  while (start[end]>='0' && start[end]<='9') {
    x = x*10+(start[end]-'0');
    ++end;
  }
  ++end; /* skip '|' */
  return end - old;
}

/*return the read bytes. */
long getLong(long& offset, long& p1){
  long x = 0;
  long old = offset;
  while (start[offset]>='0' && start[offset]<='9') {
    x = x*10+(start[offset]-'0');
    ++offset;
  }
  ++offset; /* skip '\n' */
  p1 = x;
  return offset-old;
}

void getString(char* name){
    int len = 0;
    while(start[end] != '|'){
        name[len] = start[end];
        ++len;
        ++end;
    }
    name[len] = 0;
    ++end; /* skip '|' */
}
};
