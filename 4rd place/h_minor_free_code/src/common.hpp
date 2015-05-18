#ifndef GUARD_COMMON
#define GUARD_COMMON

#include <string>
#include <iostream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include "osv.hpp"

typedef long long LL;
typedef unsigned long long ULL;

inline long long get_fsize(string fname) {
  int fd = open(fname.c_str(),O_RDONLY);
	struct stat sb;
	int err;
	err = fstat(fd, &sb);
	if (err) {
		fprintf(stderr, "fstat error! [%s]\n", strerror(errno));
		return -1;
	}
  close(fd);
	return sb.st_size;
}

double get_current_time_sec(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

void message(int query_type, const std::string &msg, double start, bool debug = true){
  if (debug){
    cerr << "Q" << query_type << ")"
         << std::setw(6)
         <<(long)((get_current_time_sec() - start)*1000) << "ms "
         << msg << std::endl;
  }
}

void load_pairs(const std::string &file, std::vector<std::pair<LL, LL> > &ps,int cap = 0){
  std::vector<std::pair<LL, LL> > my_ps;
  my_ps.reserve(cap);
  OrSeparatedValues osv(file);
  osv.nextline();

  const char * buf = osv.mmapbuff + osv.bidx;
  const char * buf_end = osv.mmapbuff + osv.fsize;
  while (buf < buf_end){
    char ch;
    LL ret;
    ret = (*buf++) & 15;
    while((ch = (*buf++)) != '|') ret = ret * 10 + (ch & 15);
    LL u = ret;
    ret = (*buf++) & 15;
    while((ch = (*buf++)) != '\n') ret = ret * 10 + (ch & 15);
    LL v = ret;
    my_ps.push_back(make_pair(u, v));
  }
  my_ps.swap(ps);
}

#endif
