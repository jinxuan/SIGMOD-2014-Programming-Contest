#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cctype>
#include <cassert>
#include <algorithm>
#include <queue>
#include <map>
#include <sys/time.h>
#include <unordered_map>
#include <tuple>
#include "union_find.hpp"
#include "osv.hpp"
#include "common.hpp"
using namespace std;

const bool debug1 = true;

struct edge_t{
  int person1;
  int person2;
  int freq;
  int rev;
  edge_t(){}
  edge_t(int person1, int person2, int freq, int rev = -1)
    : person1(person1), person2(person2), freq(freq), rev(rev) {}
};

struct less_index{
  bool operator()(const edge_t &e1, const edge_t &e2){
    if (e1.person1 == e2.person1){
      return e1.person2 < e2.person2;
    } else {
      return e1.person1 < e2.person1;
    }
  }
};

struct more_frequency{
  bool operator()(const edge_t &e1, const edge_t &e2){
    return e1.freq > e2.freq;
  }
};

struct query_t{
  int person1;
  int person2;
  int threshold;
  int index;
  query_t(int person1, int person2, int threshold, int index) :
    person1(person1), person2(person2), threshold(threshold), index(index){}
  
  bool friend operator<(const query_t &q1, const query_t &q2){
    return q1.threshold < q2.threshold;
  }
  bool friend operator>(const query_t &q1, const query_t &q2){
    return q1.threshold > q2.threshold;
  }
};

int read_number(const string &str, size_t &cur){
  int sign   = 1;
  int number = 0;
  while (str[cur] == ' ') cur++;
  if (str[cur] == '-'){
    sign = -1;
    cur++;
  }
  while (isdigit(str[cur])) number = number * 10 + str[cur++] - '0';
  return number * sign;
}

int load_edges(const string &graph_file, vector<edge_t> &es){

  int n = 0;
  vector<edge_t> my_es;
  OrSeparatedValues osv(graph_file);
  osv.nextline();
  
  const char * buf = osv.mmapbuff + osv.bidx;
  const char * buf_end = osv.mmapbuff + osv.fsize;
  while (buf < buf_end){
    char ch;
    int ret;
    ret = (*buf++) & 15;
    while((ch = (*buf++)) != '|') ret = ret * 10 + (ch & 15);
    int u = ret;
    ret = (*buf++) & 15;
    while((ch = (*buf++)) != '\n') ret = ret * 10 + (ch & 15);
    int v = ret;
    if(v > n) n = v;
    if(u < v) my_es.push_back(edge_t(u, v, 0));
  }
  n++;
  sort(my_es.begin(), my_es.end(), less_index());
  my_es.swap(es);
  return n;
}


class PairsLoader{
public:
  static const long long mask=0x0f0f0f0f0f0f0f0fLL;
  OrSeparatedValues osv;
  const char *buf;
  const char *buf_end;
  int digit_p1;

  PairsLoader(const string& fileName)
    :osv(fileName),digit_p1(1){
    
    osv.nextline();
    buf = osv.mmapbuff + osv.bidx;
    buf_end = osv.mmapbuff + osv.fsize;
  }
  
  ~PairsLoader(){
  }
  
  inline long long loadulong(char sep){
    if('0'<=buf[0] && buf[0]<='9'){
    }else{cerr<<(int)buf[0]<<"bubububub\n"<<endl;
    }
    char ch;
    long long ret = (*buf++) & 15;
    while((ch = (*buf++)) != sep) ret = ret * 10 + (ch & 15);
    return ret;
  }
  
  inline bool is_end(){
    return buf >= buf_end;
  }
  
  __attribute__((always_inline))  inline pair<long long,int> loadfastnum(char sep){
    assert(!is_end());
    if(buf + digit_p1 >= buf_end)goto init_digit;
    if(buf[digit_p1]!=sep){
      digit_p1++;
      if(buf + digit_p1 >= buf_end || buf[digit_p1]!=sep){
      init_digit:;
        digit_p1=1;
        while(buf[digit_p1] != sep)digit_p1++;
      }
      assert(digit_p1 <= (15+1));

    }
    
    long long ret;
    int lower;
    if(digit_p1>=9){
      int shift = (8+8)-digit_p1;
      unsigned long long num_u = __builtin_bswap64(((const unsigned long long *)(buf-shift))[0])&mask;
      unsigned long long num_l = __builtin_bswap64(((const unsigned long long *)(buf-shift))[1])&mask;
      num_u&=(~ 0ULL)>>(shift*8);
      num_u<<=4;
      ret = num_u|num_l;
      lower = num_l;
    }else{
      int shift = 8-digit_p1;
      unsigned long long num = __builtin_bswap64(((const unsigned long long *)(buf-shift))[0]);
      ret = num&mask&((~ 0ULL)>>(shift*8));
      lower = ret;
    }
    buf+=digit_p1+1;
    lower |= lower>>4;
    lower &= 0xff00ff;
    return make_pair(ret,((lower>>8)|lower)&0xffff);
  }
};


void load_graph_fast(const string &graph_dir, vector<edge_t> &es){
  string tmp;
  double s = get_current_time_sec();
  int    n = 0;
  vector<vector<int> > cr_vec_tmp_g;
  
  
  s = get_current_time_sec();

  {
    FILE* fp=fopen("endq1load","w");
    if(fp) fclose(fp);
  }
  
#pragma omp parallel sections
  {
#pragma omp section
    {
      vector<vector<int> > cr_vec_tmp(10000);
      string comment_creator = graph_dir + "/comment_hasCreator_person.csv";
      string comment_reply = graph_dir + "/comment_replyOf_comment.csv";
      const int buff_size = 1 << 20;
      const int buff_mask = buff_size - 1;
      int buff[buff_size];

      memset(buff, 0, sizeof(buff));
      PairsLoader cpl(comment_creator);
      PairsLoader crl(comment_reply);
      LL last_r = -1;
      int p_id  = -1;
      
      while(!crl.is_end()){
        const auto r_raw = crl.loadfastnum('|');
        const auto c_raw = crl.loadfastnum('\n');
        
        LL r = r_raw.first;
        while (last_r != r){
          auto cp_raw = cpl.loadfastnum('|');
          p_id        = cpl.loadulong('\n');
          buff[cp_raw.second & buff_mask] = p_id;
          last_r = cp_raw.first;
        }
        
        int p1 = p_id;

        int p2 = buff[c_raw.second & buff_mask];

        if(max(p1,p2) >= (int)cr_vec_tmp.size()){
          int num = 10;
          while(max(p1,p2) >= (int)cr_vec_tmp.size()*num)num*=10;
          vector<vector<int> > cr_vec_tmp_new(cr_vec_tmp.size() * num);
          for(int i = 0; i < (int)cr_vec_tmp.size(); i++) cr_vec_tmp_new[i].swap(cr_vec_tmp[i]);
          cr_vec_tmp_new.swap(cr_vec_tmp);
        }
        if(p1 < p2){
          cr_vec_tmp[p1].push_back(p2*2);
        }else{
          cr_vec_tmp[p2].push_back(p1*2+1);
        }
      }
      cr_vec_tmp_g.swap(cr_vec_tmp);
    }
#pragma omp section
    {
      n = load_edges(graph_dir + "/person_knows_person.csv", es);
    }
  }
  
  message(1, "comment->person", s, debug1);
  
  s = get_current_time_sec();

  {
    vector<tuple<int,int,int>> cr_vec(n);
    size_t es_i = 0;
    for (int i = 0; i < n; i++){
      for (auto it = cr_vec_tmp_g[i].begin(); it != cr_vec_tmp_g[i].end(); it++){
        int itn = *it;
        auto & v= cr_vec[itn>>1];
        
        if(get<0>(v) != i) v = make_tuple(i,0,0);//initialize

        if(itn&1){
          get<1>(v)++;
        } else {
          get<2>(v)++;
        }
      }
      for (; es_i < es.size() && es[es_i].person1 == i; es_i++){
        const auto &v = cr_vec[es[es_i].person2];
        int freq = get<0>(v) != i ? 0 : min(get<1>(v),get<2>(v));
        es[es_i].freq = freq;
      }
    }
  }
  message(1, "person's comment & Construct graph",  s, debug1);
}

void load_query(const string &query_file, vector<query_t> &qs){
  
  string   str;
  ifstream ifs(query_file.c_str());
  
  if (!ifs.is_open()){
    perror("ifstream");
    exit(EXIT_FAILURE);
  }

  while (getline(ifs, str)){
    size_t cur = 0;
    for (size_t i = 0; i < str.size(); i++)
      if (!isdigit(str[i]) && str[i] != '-') str[i] = ' ';
    
    if (read_number(str, cur) == 1){
      int p1    = read_number(str, cur);
      int p2    = read_number(str, cur);
      int thres = read_number(str, cur);
      qs.push_back(query_t(p1, p2, thres, qs.size()));
    }
  }
  
  ifs.close();
}

int distance(const vector<vector<int> > &adj_list,
             int p1,
             int p2,
             vector<bool> &visit1,
             vector<bool> &visit2){
  int res = -1;
  int d   = 0;
  int c   = 0;
  
  if (p1 == p2) return d;
  queue<int> que1[2];
  queue<int> que2[2];
  vector<int> updated1, updated2;

  que1[0].push(p1);
  que2[0].push(p2);
  updated1.push_back(p1);
  updated2.push_back(p2);
  visit1[p1] = true;
  visit2[p2] = true;
  
  while (!que1[c].empty() || !que2[c].empty()){
    d++;
    while (!que1[c].empty()){
      int v = que1[c].front(); que1[c].pop();
      for (size_t i = 0; i < adj_list[v].size(); i++){
        int w = adj_list[v][i];
        if (visit2[w]){
          res = d;
          goto END;
        }
        if (!visit1[w]){
          visit1[w] = true;
          que1[1 - c].push(w);
          updated1.push_back(w);
        }
      }
    }
    d++;
    while (!que2[c].empty()){
      int v = que2[c].front(); que2[c].pop();
      for (size_t i = 0; i < adj_list[v].size(); i++){
        int w = adj_list[v][i];
        if (visit1[w]){
          res = d;
          goto END;
        }
        if (!visit2[w]){
          visit2[w] = true;
          que2[1 - c].push(w);
          updated2.push_back(w);
        }
      }
    }
    c = 1 - c;
  }
 END:
  for (size_t i = 0; i < updated1.size(); i++) visit1[updated1[i]] = false;
  for (size_t i = 0; i < updated2.size(); i++) visit2[updated2[i]] = false;
  return res;
}
             
int main(int argc, char **argv){
  double start = get_current_time_sec();
  ios::sync_with_stdio(false);
  
  if (argc != 3){
    cerr << "Usage: " << argv[0]
         << " (path_to_data_directory)"
         << " (path_to_query_file)" << endl;
    exit(EXIT_FAILURE);
  }
  
  vector<edge_t>  es;
  vector<query_t> qs;
  load_graph_fast(argv[1], es);
  load_query(argv[2], qs);
  sort(es.begin(), es.end(), more_frequency());
  sort(qs.begin(), qs.end(), greater<query_t>());
  

  int V = 0;
  for (size_t i = 0; i < es.size(); i++){
    V = max(V, es[i].person1 + 1);
    V = max(V, es[i].person2 + 1);
  }
  
  vector<vector<int> > adj_list(V);

  size_t pos = 0;
  vector<int > q_ans(qs.size(), -1);
  vector<bool> visit1(V, false);
  vector<bool> visit2(V, false);
  UnionFind    uf(V);
  
  
  {
    double s = get_current_time_sec();
    
    for (size_t i = 0; i < qs.size(); i++){
      while (pos < es.size() && es[pos].freq > qs[i].threshold){
        adj_list[es[pos].person1].push_back(es[pos].person2);
        adj_list[es[pos].person2].push_back(es[pos].person1);
        uf.unite(es[pos].person1, es[pos].person2);
        pos++;
      }
      if (uf.same(qs[i].person1, qs[i].person2)){
        q_ans[qs[i].index] =
          distance(adj_list, qs[i].person1, qs[i].person2, visit1, visit2);
      }
    }
    message(1, "Execute queries", s, debug1);
  }
  
  for (size_t i = 0; i < q_ans.size(); i++){
    cout << q_ans[i] << endl;
  }
  message(1, "Total", start, true);
  exit(EXIT_SUCCESS);
}

