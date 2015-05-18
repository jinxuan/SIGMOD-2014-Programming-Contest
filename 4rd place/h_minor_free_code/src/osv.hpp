#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <smmintrin.h>
#include <assert.h>
#include <xmmintrin.h>

#define bsf(x) __builtin_ctzll(x)

using namespace std;

inline long long get_fsize(int fd) {
  struct stat sb;
  int err;
  err = fstat(fd, &sb);
  if (err) {
    fprintf(stderr, "fstat error! [%s]\n", strerror(errno));
    return -1;
  }
  return sb.st_size;
}

class OrSeparatedValues{
  FILE* fp;
  int fd;
public:
  const char* mmapbuff;
  long long bidx;
  long long csize;
  long long fsize;

  inline void ungetchar(){
    --bidx;
  }

  OrSeparatedValues(const string& fileName):
    bidx(0),csize(0)
  {
    fd = open(fileName.c_str(), O_RDONLY);
    if(fd == -1){
      cerr << "Open file " << fileName << " failed" << endl;
      exit(EXIT_FAILURE);
    }
        
    fsize =  get_fsize(fd);

    if(fsize < 0){
      cerr << "getfilesize error" << endl;
      exit(EXIT_FAILURE);
    }

    mmapbuff = (const char*)mmap(NULL, fsize, PROT_READ, MAP_SHARED, fd, 0);
  }

  ~OrSeparatedValues(){
    int err = munmap((void *)mmapbuff, fsize);
    if (err) {
      fprintf(stderr, "munmap error! [%d]\n", err);
    }
    close(fd);
  }

  int my_getchar(){
    if(bidx == fsize) return -1;
    int ch = mmapbuff[bidx++];
    return ch;
  }
  
  int my_getchar_simple(){
    return mmapbuff[bidx++];
  }
  
  bool nextline(){
    while(true){
      int ch = my_getchar();
      if(ch == -1) return false;
      if(ch == '\n') break;
      if(ch == '\r'){
        ch = my_getchar();
        if(ch != '\n' && ch!=-1) ungetchar();
        break;
      }
    }
    int ch = my_getchar();
    if(ch == -1) return false;
    ungetchar();
    return true;
  }

  //assume newline as \r\n or \n
  bool nextline_maybefar(){
    const __m128i allcharLF = _mm_set1_epi8('\n');
    while(true){
      if(fsize -bidx < 48)return nextline();
      const char * current = mmapbuff +bidx;
      const char * base = (char *)(((size_t)current) & ~((size_t)15));
      int diff = current - base;
		  
      __m128i mask1 = _mm_cmpeq_epi8(allcharLF,_mm_load_si128((const __m128i *)(base)));
      __m128i mask2 = _mm_cmpeq_epi8(allcharLF,_mm_load_si128((const __m128i *)(base+16)));
      __m128i mask3 = _mm_cmpeq_epi8(allcharLF,_mm_load_si128((const __m128i *)(base+32)));
      unsigned long long maskbit1 = _mm_movemask_epi8(mask1);
      unsigned long long maskbit2 = _mm_movemask_epi8(mask2);
      unsigned long long maskbit3 = _mm_movemask_epi8(mask3);
		  
      unsigned long long all = ((maskbit3<<32)|(maskbit2<<16)|maskbit1)>>diff;
		  
      int pos = bsf(all);
		  
      if(all==0){
        bidx+=48-diff;
        continue;
      }else{
        bidx+=pos+1;
      }
      return bidx < fsize;
    }
  }
	
  int nextuint(){
    char ch;
    while(!isdigit(ch = my_getchar()));
    int ret=ch&15;
    while(isdigit(ch = my_getchar()))
      ret = (ret<<3)+(ret<<1)+(ch&15);
    if(ch != -1) ungetchar();
    return ret;
  }

  long long nextulong(){
    char ch;
    while(!isdigit(ch = my_getchar()));
    long long ret=ch&15;
    while(isdigit(ch = my_getchar()))
      ret = (ret<<3)+(ret<<1)+(ch&15);
    if(ch != -1) ungetchar();
    return ret;
  }
  
  long long nextulong1(){
    char ch;
    
    long long ret= my_getchar_simple() & 15;
    
    while((ch = my_getchar_simple()) != '|')
      ret = ret * 10 + (ch & 15);

    return ret;
  }
  
  bool is_end(){
    return bidx >= fsize;
  }

  void next_bar(){
    while (my_getchar_simple() != '|');
  }

  bool q4_nextline(){
    bidx += 21;
    return bidx < fsize;
  }

  string next(){
    string ret;
    int ch = my_getchar();
    if(ch == -1) return ret;
    if(ch != '|' && ch != '\n' && ch != '\r') ret = ch;

    while(true){
      ch = my_getchar();
      if(ch == -1 || ch == '|' || ch == '\n' || ch == '\r') return ret;
      ret += (char)ch;
    }
  }
};


class OrSeparatedValuesPairs{
  static inline long long to_num(__m128i num)
  {
    const __m128i num10 = _mm_set1_epi16(10 * 256+1);
    const __m128i num100 = _mm_set1_epi32(100 * 256 * 256+1);
    const __m128i num10000 = _mm_set1_epi64x(10000LL * 256 * 256 * 256 * 256+1);

    __m128i num2 = _mm_maddubs_epi16(num, num10);
    __m128i num3 = _mm_madd_epi16(num2, num100);
    __m128i num4 = _mm_mullo_epi32(num3, num10000);
    __m128i num5 = _mm_hadd_epi32(num4, num4);
    return _mm_extract_epi32(num5, 1) * (long long)100000000 + _mm_extract_epi32(num5, 0);
  }
  static inline __m128i make_mask(__m128i mask)
  {
    mask = _mm_or_si128(mask, _mm_slli_si128(mask, 1));
    mask = _mm_or_si128(mask, _mm_slli_si128(mask, 2));
    mask = _mm_or_si128(mask, _mm_slli_si128(mask, 4));
    return _mm_or_si128(mask, _mm_slli_si128(mask, 8));
  }
  int fd;
  const char* mmapbuff;
  long long fsize;

public:
  OrSeparatedValuesPairs(const string& fileName)
  {
    fd = open(fileName.c_str(), O_RDONLY);
    if(fd == -1){
      cerr << "Open file " << fileName << " failed" << endl;
      exit(EXIT_FAILURE);
    }

    fsize =  get_fsize(fd);

    if(fsize < 0){
      cerr << "getfilesize error" << endl;
      exit(EXIT_FAILURE);
    }

    mmapbuff = (const char*)mmap(NULL, fsize, PROT_READ, MAP_SHARED, fd, 0);
  }

  ~OrSeparatedValuesPairs(){
    int err = munmap((void *)mmapbuff, fsize);
    if (err) {
      fprintf(stderr, "munmap error! [%d]\n", err);
    }
    close(fd);
  }
  
  void load_pairs(vector<pair<int, int> > &pairs,bool need_reverse=true)
  {
    const char* mmapbuff = this->mmapbuff;
    long long fsize = this->fsize;

    long long sidx = 0;
    while (sidx < fsize && mmapbuff[sidx] != '\n')sidx++;
    while (sidx < fsize && !isdigit(mmapbuff[sidx]))sidx++; 

    if (sidx >= fsize)return; //empty file
    assert(isdigit(mmapbuff[sidx]));
    sidx -= 15;
    assert(sidx >= 0);

    const __m128i allchar0 = _mm_set1_epi8('0');
    const __m128i allchar0_m = _mm_set1_epi8('0'+128);
    const __m128i allchar9_m = _mm_set1_epi8(-128+9);
    const __m128i reverse_arr = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);

    long long idx = fsize - 16 - 1;


    while (idx > sidx) { 
      __m128i val1_raw,val2_raw;
      {
      load_val1:
        __m128i buf = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)(mmapbuff + idx)), reverse_arr);
        __m128i num = _mm_sub_epi8(buf, allchar0);
        __m128i mask = _mm_cmpgt_epi8(_mm_sub_epi8(buf, allchar0_m), allchar9_m);
        unsigned short maskbit = _mm_movemask_epi8(mask);
        assert(maskbit !=0);
        if (maskbit & 1) {
          idx -= bsf(maskbit + 1);
          goto load_val1;
        }else {
          mask = make_mask(mask);
          val1_raw = _mm_andnot_si128(mask, num);
          idx -= bsf(maskbit) + 1;
        }
      }
      {
        __m128i buf = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)(mmapbuff + idx)), reverse_arr);
        __m128i num = _mm_sub_epi8(buf, allchar0);
        __m128i mask = _mm_cmpgt_epi8(_mm_sub_epi8(buf, allchar0_m), allchar9_m); 
        unsigned short maskbit = _mm_movemask_epi8(mask);
        assert(maskbit!=0);
        mask = make_mask(mask);
        val2_raw = _mm_andnot_si128(mask, num);
        idx -= bsf(maskbit) + 1;
      }
      int val1 = to_num(val1_raw);
      int val2 = to_num(val2_raw);
      pairs.push_back(make_pair(val2,val1));
    }
    if(need_reverse) reverse(pairs.begin(),pairs.end());
  }
};
