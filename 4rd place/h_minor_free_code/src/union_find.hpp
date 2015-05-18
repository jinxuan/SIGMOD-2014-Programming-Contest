/*********************************************************
UnionFind

find  : O(logn)
weight: O(logn)
same  : O(logn)
unite : O(logn)

verified at
http://www.spoj.com/problems/MST/ (only find, same and unite are used.)
***********************************************************/

#ifndef GUARD_UNION_FIND
#define GUARD_UNION_FIND


#include <vector>

class UnionFind{
  std::vector<int> parent;
  std::vector<int> count;
  std::vector<int> rank;
public:
  UnionFind(int N) : parent(std::vector<int>(N)),
                     count(std::vector<int>(N, 1)),
                     rank(std::vector<int>(N, 0)){
    for(int i = 0; i < N; i++) parent[i] = i;
  }
  
  int find(int x){
    if(x == parent[x]) return x;
    else return parent[x] = find(parent[x]);
  }
  
  int weight(int x){
    return count[find(x)];
  }
  
  bool same(int x, int y){
    return find(x) == find(y);
  }
    
  bool unite(int x, int y){
    x = find(x);
    y = find(y);
    if(x == y) return false;
    if(rank[x] < rank[y]){
      count[y] += count[x];
      parent[x] = y;
    }else{
      count[x] += count[y];
      parent[y] = x;
      if(rank[x] == rank[y]) rank[y]++;
    }
    return true;
  }
};


#endif
