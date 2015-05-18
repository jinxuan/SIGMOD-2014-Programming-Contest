#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
#include "common.hpp"
#include "union_find.hpp"
#include "osv.hpp"
#include <omp.h>

using namespace std;
typedef long long LL;

namespace {
  const bool debug4  = true;
};


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

bool exist_file(string a){
  struct stat t;
  return !stat(a.c_str(),&t);
}

/* 
   To speed up bfs, we iterate edges from unreached vertices
   after visiting some fraction of vertices from a source vertex.
*/
int backward_bfs(const vector<int> &vs,
                 const int *edge,
                 const vector<int> &edge_start,
                 vector<char> visit,
                 int reached_dist,
                 int *array_que,
                 int dsum,
                 int max_dsum)
{
  int  cc    = reached_dist;
  int  cur   = 0;
  int *cur_p = array_que;
  int *end_p = array_que;
  
  for (auto it = vs.begin(); it != vs.end(); it++){
    *end_p++ = *it;
    end_p -= visit[*it];
  }
  dsum += cc*(end_p-array_que);
  
  while (cur_p != end_p){
    cc++;
    cur++;
    dsum += (end_p-array_que);
    int *nxt_p = array_que;
    for (const int *p = cur_p; p != end_p; p++){
      int  cv = *p;
      bool ok = false;

      const int *edge_ptr = edge + edge_start[cv];
      const int *edge_end = edge + edge_start[cv + 1];
      
      while (!ok && edge_ptr != edge_end){
        ok = ok || visit[*edge_ptr++] == cur;
      }
      
      if (ok){
        visit[cv] = cur + 1;
      } else {
        *nxt_p++ = cv;
        int num_next    = nxt_p - array_que;
        if(dsum + num_next > max_dsum)return 0;
      }
    }
    end_p = nxt_p;
  }
  return dsum;
}

/*
  compute the lowerbound of sum of distance from a source vertex s by estimating
  the number of vertices whose distance from s is not larger than dist + 1
*/

inline int estimate_lowerbound(const int *que_begin,
                               const int *que_end,
                               int dist,
                               int current_distance_sum,
                               int num_of_unreached_nodes,
                               const vector<int> &edge_start){
  int approx_next1 = 0;
  for (const int *p = que_begin; p != que_end; p++){
    approx_next1 += edge_start[(*p) + 1] - edge_start[*p] - 1;
  }
  
  if (approx_next1 <= num_of_unreached_nodes){
    int plus_one = approx_next1;
    int plus_two = num_of_unreached_nodes - approx_next1;
    return current_distance_sum + plus_one * (dist + 1) + plus_two * (dist + 2);
  } else {
    return -1e9;
  }
}

/*
  compute the top-k central people exactly by conducting bfs with pruning.
  openmp is used to speed up computation
*/
void pruned_bfss(int n,
                 int k,
                 int component_size,
                 const vector<vector<int> > &es,
                 const vector<int> &vs,
                 set<pair<int, int> > &answer,
                 int max_dsum){
  
  assert(component_size == (int)vs.size());
  
  // assume no self-loop in graph
  const int inf      = 1e9;
  int edge_num = 0;

  for (size_t i = 0; i < es.size(); i++){
    edge_num += es[i].size();
  }

  vector<int> edge_vec;
  vector<int> edge_start;
  edge_vec.reserve(edge_num);
  edge_start.reserve(es.size() + 1);
  edge_start.push_back(0);
  
  {
    int start_index = edge_start[0] = 0;
    for (size_t i = 0; i < es.size(); i++){
      edge_vec.insert(edge_vec.end(),es[i].begin(),es[i].end());
      start_index += es[i].size();
      edge_start.push_back(start_index);
    }
  }
  const int *edge = &(edge_vec[0]);

#ifdef _OPENMP
  int thread_num = 6;
  if(exist_file("endq1")){
    cerr<<"e";
    ++thread_num;
  }
  
  if(exist_file("endq1load")){
    cerr<<"l";
    ++thread_num;
  }
  cerr<<thread_num<<"t";
  
  omp_set_num_threads(thread_num);  
#endif
  
  vector<int> dsum_arr(es.size());
  
  
  for (int ss = 0; ss < (int)vs.size();){
    int limit = min(max(ss * 8 + 7, k), (int)vs.size() - ss);
    if(ss > (int)vs.size()/100)
      limit = vs.size() - ss;
    vector<pair<int, int> > vec_ans(limit, make_pair(inf, 0));
    
#pragma omp parallel 
    {
      vector<char> visit(n);
      vector<int> array_que_vec(component_size);
      int *array_que = &(array_que_vec[0]);
      
#pragma omp for nowait schedule(dynamic)

      for(int para = 0; para < limit; para++){
        const int s    = vs[ss + para];

        if(edge_start[s+1] - edge_start[s] == 1){
          int next = edge[edge_start[s]];
          if(dsum_arr[next] > max_dsum) continue;
        }

        int dsum = 0;
      
        visit.assign(n, 0);
      
        int* cur_p = array_que;
        int* nxt_p = array_que;
        int* end_p = array_que;
      
        for (const int *p = edge + edge_start[s]; p != edge + edge_start[s + 1]; p++){
          *end_p++  = *p;
          visit[*p] = 1;
        }
      
        visit[s] = 1;
        int num_of_unreached_nodes = component_size - 1;
        LL  cc                     = 1;
      
        while(nxt_p != end_p){
          cur_p = nxt_p;
          nxt_p = end_p;
        
          num_of_unreached_nodes -= nxt_p - cur_p;
          dsum                   += cc * (nxt_p - cur_p);
        
          if(dsum + num_of_unreached_nodes * (cc + 1) > max_dsum) goto CONT;
        
          int dsum_lb = estimate_lowerbound(cur_p, end_p, cc, dsum,
                                            num_of_unreached_nodes, edge_start);
          if (dsum_lb > max_dsum) goto CONT;
        
          if ((nxt_p - cur_p)*10 > num_of_unreached_nodes){
            dsum = backward_bfs(vs, edge, edge_start, visit, cc, array_que, dsum, max_dsum);
            if(dsum == 0)goto CONT;
            goto END;
          }
        
          while(cur_p != nxt_p){
            int cv = *cur_p++;

            const int *edge_ptr = edge + edge_start[cv];
            const int *edge_end = edge + edge_start[cv + 1];
            while (edge_ptr != edge_end){
              const int to = *edge_ptr++;
              *end_p = to;
              end_p  += 1 - visit[to];
              visit[to] = 1;
            }
            dsum_lb += edge_start[cv+1]- edge_start[cv] - 1;
          
            if (dsum_lb - (end_p - nxt_p) > max_dsum){
              goto CONT;
            } else if(num_of_unreached_nodes == end_p - nxt_p){
              dsum += (cc + 1) * (end_p - nxt_p);
              goto END;
            }
          }
          ++cc;
        }
      END:;
        vec_ans[para] = make_pair(dsum, s);

        dsum_arr[s] = dsum;
        continue;
      CONT:;
        dsum_arr[s] = inf;
      }
    }
    
    for(int para = 0; para < limit; para++){
      answer.insert(vec_ans[para]);
      if((int)answer.size() > k){
        answer.erase(*answer.rbegin());
      }
      if((int)answer.size() >= k){
        max_dsum = min(max_dsum ,answer.rbegin()->first);
      }
    }
    ss += limit;
  }
}

// compute top-k central people for each connected component.
vector<pair<double, int> > compute_top_k(int n,
                                         int k,
                                         const vector<vector<int> > &es,
                                         const vector<int> &vs,
                                         double max_cent){
  const int inf      = 1e9;
  int       max_dsum = max_cent==0 ? inf : min(ceil((double)(vs.size() - 1) * (vs.size() - 1) / (n - 1) / max_cent),1e9);
                                          

  set<pair<int, int> > answer;
  pruned_bfss(n, k, vs.size(), es, vs, answer, max_dsum);
  
  vector<pair<double, int> > top_k;
  for (auto it = answer.begin(); it != answer.end(); it++){
    int    dsum = it->first;
    double cent = (double)(vs.size() - 1) * (vs.size() - 1) / (n - 1) / dsum;
    top_k.push_back(make_pair(-cent, it->second));
  }
  return top_k;
}

void query(const vector<int> &_vs,
           const vector<vector<int> > &adj,
           int k,
           string& out,
           vector<int> &vertex_id) {
    
  int n = _vs.size();
  vector<vector<int> > es(n);
  UnionFind uf(n);
  
  if(!_vs.empty()){
    
    for (size_t i = 0; i < _vs.size(); i++){
      vertex_id[_vs[i]] = i;
    }

    int count = _vs.size() - 1;

    for (int i = 0; i < (int)_vs.size(); i++){
      int v = _vs[i];
      
      for (auto it = adj[v].begin(); it != adj[v].end(); it++){
        int j = vertex_id[*it];
        if (j == -1) continue;
        es[i].push_back(j);
        if (i < j && count > 0){
          if (uf.unite(i, j)) count--;
        }
      }
    }
    
    for (size_t i = 0; i < _vs.size(); i++){
      vertex_id[_vs[i]] = -1;
    }
  
  }
  
  set<pair<double, int> > answer;
  vector<pair<int, int> > sorted_by_degree;

      
  for(int s = 0; s < n; ++s){
    sorted_by_degree.push_back(make_pair(-(int)es[s].size(), s));
  }
  
  sort(sorted_by_degree.begin(), sorted_by_degree.end());
  
  vector<vector<int> > connected_components(n);

  for (size_t i = 0; i < sorted_by_degree.size(); i++){
    int v = sorted_by_degree[i].second;
    connected_components[uf.find(v)].push_back(v);
  }

  vector<pair<double, int> > top_k;
  
  for (int i = 0; i < n; i++){
    if ((int)connected_components[i].size() > 1){
      double max_cent = (int)top_k.size() < k ? 0 : -top_k.rbegin()->first;
      vector<pair<double, int> > component_top_k =
        compute_top_k(n, k, es, connected_components[i],max_cent);
      top_k.insert(top_k.end(), component_top_k.begin(), component_top_k.end());
      if((int)top_k.size() >= k){
        sort(top_k.begin(), top_k.end());
        top_k.erase(top_k.begin()+k,top_k.end());
      }
    }
  }

  
  for (int i = 0; i < n; i++)
    if ((int)connected_components[i].size() == 1)
      top_k.push_back(make_pair(0, connected_components[i][0]));
  sort(top_k.begin(), top_k.end());
  if ((int) top_k.size() > k) top_k.resize(k);  
  

  stringstream ss;
  if ((int)top_k.size() < k){
    cerr << "line" << __LINE__ << ": top_k.size() < k" << endl;
  }
  
  for (int i = 0; i < (int)top_k.size(); i++) {
    int x = top_k[i].second;
    if(i){
      ss << " ";
    }
    ss << _vs[x];
  }
  out = ss.str();
}

int main(int argc, char *argv[]) {

  if (argc != 3){
    cerr << "Usage: " << argv[0]
         << " (path_to_data_directory)"
         << " (path_to_query_file)" << endl;
    exit(EXIT_FAILURE);
  }

  const string data_dir   = string(argv[1]);
  const string query_file = string(argv[2]);
  double start = get_current_time_sec();

  int num_of_persons = 0;

  /* forum_hasTag_tag */
  vector<pair<LL, LL> > tag2forum;

  vector<vector<int> > adj;

  /* forum_hasMember_person */
  vector<pair<LL, LL> > forum2person;

  // read query
  vector<pair<int, LL> > queries;
  
#pragma omp parallel sections
  {
#pragma omp section
    {
      /* tag */
      map<string, LL> tagname2tagid;
      
      {
        OrSeparatedValues osv(data_dir + "tag.csv");
        for (int linum = 0; osv.nextline(); linum++) {
          LL     id   = osv.nextuint();
          string name = osv.next();
          tagname2tagid[name] = id;
        }
      }
      {
        OrSeparatedValues osv(data_dir + "forum_hasTag_tag.csv");
        for (int linum = 0; osv.nextline(); linum++) {
          LL forum = osv.nextulong();
          LL tag   = osv.nextuint();
          tag2forum.push_back(make_pair(tag, forum));
        }
        sort(tag2forum.begin(), tag2forum.end());
      }
      {
        ifstream ifs(query_file.c_str());
        char cs[1 << 10];
        
        for (string s; getline(ifs, s);) {
          if (s.find("query4(", 0) != string::npos) {
            int k;
            sscanf(s.c_str(), "query4(%d, %s)", &k, cs);
            string tag_name(cs);
            tag_name = tag_name.substr(0, tag_name.length() - 1);
            LL tag = tagname2tagid[tag_name];
            queries.push_back(make_pair(k, tag));
          }
        }
      }
      {
        const int buf_size = (1 << 20);
        const int buf_mask = buf_size - 1;
        vector<char> maybe_important_forum(buf_size, false);

        
        for (size_t q = 0; q < queries.size(); q++) {
          LL tag = queries[q].second;
          auto it = lower_bound(tag2forum.begin(), tag2forum.end(),
                                make_pair(tag, 0LL));
          for (; it->first == tag; it++) {
            maybe_important_forum[(it->second>>4) & buf_mask] = true;
          }
        }

        {
          OrSeparatedValues osv(data_dir + "forum_hasMember_person.csv");
          while(osv.nextline_maybefar()){
            LL forum  = osv.nextulong1();
            if (maybe_important_forum[(forum>>4) & buf_mask]){
              LL person = osv.nextulong1();
              forum2person.push_back(make_pair(forum, person));
            }
          }
        }
        sort(forum2person.begin(), forum2person.end());
      }
    }
    
#pragma omp section
    {
      /* person_knows_person */
      vector<pair<int,int> > links;
      OrSeparatedValuesPairs osv(data_dir + "person_knows_person.csv");
      osv.load_pairs(links,false);
      
      sort(links.begin(),links.end());
      num_of_persons = links.rbegin()->first + 1;
      
      adj.assign(num_of_persons,vector<int>());
      for (auto it = links.begin(); it != links.end(); it++){
        adj[it->first].push_back(it->second);
      }
    }
  }

  message(4, "Load files", start, debug4);
  start = get_current_time_sec();
  vector<string> query_answer(queries.size());
  
  vector<int> vertex_id(num_of_persons, -1);


  for (size_t q = 0; q < queries.size(); q++) {
    vector<char> is_target(num_of_persons);
    vector<int>  targets;
    
    int k  = queries[q].first;
    LL tag = queries[q].second;
    
    auto it = lower_bound(tag2forum.begin(), tag2forum.end(),
                          make_pair(tag, 0LL));
    for (; it->first == tag; it++) {
      LL forum = it->second;
      auto jt = lower_bound(forum2person.begin(), forum2person.end(),
                            make_pair(forum, 0LL));
      for (; jt->first == forum; jt++){
        is_target[jt->second] = true;
      }
    }
    
    for (int i = 0; i < num_of_persons; i++){
      if (is_target[i]) targets.push_back(i);
    }
    
    query(targets, adj, k, query_answer[q], vertex_id);
  }

  for(size_t i = 0;i<query_answer.size(); ++i)
    cout << query_answer[i] << endl;
  message(4, "Queries", start, true);
  return 0;
}
