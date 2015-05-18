#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include "osv.hpp"
#include "common.hpp"
#include "union_find.hpp"

#define rep(i,n) for(int i=0;i<(int)(n);i++)

using namespace std;

typedef long long ll;

namespace {
  const bool debug2 = false;
  int N = 0, T = 0, Q = 0; // #vertices, #tags, #queries
  vector<int> birthday;
  vector<pair<int, int> > links;
  map<ll, int> tagToId;
  vector<string> tagName;
  vector<vector<LL> > vertexTag;
  vector<int> vertexTagSize;
  vector<pair<int, int> > queries;
  vector<vector<string> > result;
  vector<int> tagCount;
}

inline int getBirthdayCode(int y, int m, int d) {
  return y * 12 * 31 + m * 31 + d;
}

void inputFast(const string workingDir, const string queryFile);
void query();


int main(int argc, char** argv) {
  if (argc != 3) {
    exit(EXIT_FAILURE);
  }
  
  const string workingDir(argv[1]), queryFile(argv[2]);
  double t1 = get_current_time_sec();

  inputFast(workingDir, queryFile);
  message(2, "Input", t1, debug2);
  
  result = vector<vector<string> >(Q);

  double t2 = get_current_time_sec();
  query();
  message(2, "Query", t2, debug2);
  
  rep(i, Q) {
    rep(j, result[i].size()) {
      if (j) cout << " ";
      cout << result[i][j];
    }
    cout << endl;
  }
  message(2, "Total", t1, true);
  return 0;
}

void query() {
  vector<pair<pair<int, int>, int> > sortedQuery(Q);
  int kMax = 0;
  rep(i, Q) {
    sortedQuery[i] = make_pair(queries[i], i);
    kMax = max(kMax, queries[i].second);
  }
  sort(sortedQuery.begin(), sortedQuery.end(), greater<pair<pair<int, int>, int> >());
  
  vector<pair<int, pair<int, int> > > sortedEdge;
  
  for (auto it = links.begin(); it != links.end(); it++){
    int s = it->first;
    int t = it->second;
    int appearedTime = min(birthday[s], birthday[t]);
    sortedEdge.push_back(make_pair(appearedTime, make_pair(s, t)));
  }
  
  for (int i = 0; i < N; i++){
    sortedEdge.push_back(make_pair(birthday[i], make_pair(i, i)));
  }
                                                    
  sort(sortedEdge.begin(), sortedEdge.end(), greater<pair<int, pair<int, int> > >());
  

  UnionFind uf(vertexTagSize[N]);
  vector<pair<int, string> > maxComponent(T);
  vector<pair<int, string> > topKComponent(kMax);
  
  rep(i, T) {
    maxComponent[i] = make_pair(0, tagName[i]);
  }
  
  int    q = 0;
  size_t i = 0;
  while (q < Q) {
    while (i < sortedEdge.size() && sortedEdge[i].first >= sortedQuery[q].first.first) {
      int u = sortedEdge[i].second.first, v = sortedEdge[i].second.second;
      
      size_t pos_u = 0;
      size_t pos_v = 0;
      
      while (pos_u < vertexTag[u].size() && pos_v < vertexTag[v].size()){

        if (tagCount[vertexTag[u][pos_u]] < -topKComponent[kMax - 1].first){
          break;
        }
        
        if (vertexTag[u][pos_u] == vertexTag[v][pos_v]){
          uf.unite(vertexTagSize[u] + pos_u, vertexTagSize[v] + pos_v);
          
          int tag           = vertexTag[u][pos_u];
          int componentSize = uf.weight(vertexTagSize[u] + pos_u);
          
          if (componentSize > -maxComponent[tag].first) {
            maxComponent[tag].first = -componentSize;
            string componentName = tagName[tag];
            
            pair<int, string> e(-componentSize, componentName);
            
            if (e < topKComponent[kMax - 1]) {
              rep(j, kMax) {
                if (topKComponent[j].second == componentName) {
                  swap(e, topKComponent[j]);
                  break;
                }
                if (e < topKComponent[j]) {
                  swap(e, topKComponent[j]);
                }
              }
            }
          }
          pos_u++;
          pos_v++;
        } else if (vertexTag[u][pos_u] < vertexTag[v][pos_v]){
          pos_u++;
        } else {
          pos_v++;
        }
      }
      i++;
    }

    while (q < Q && (i == sortedEdge.size() || sortedEdge[i].first < sortedQuery[q].first.first)) {
      rep(j, sortedQuery[q].first.second) {
        result[sortedQuery[q].second].push_back(topKComponent[j].second);
      }
      q++;
    }
  }
}

void inputFast(const string workingDir, const string queryFile) {
  const string personFile    = workingDir + "person.csv";
  const string graphFile     = workingDir + "person_knows_person.csv";
  const string tagListFile   = workingDir + "tag.csv";
  const string vertexTagFile = workingDir + "person_hasInterest_tag.csv";
  
  {
    vector<pair<int, int> > tmpBirthday;
    
    OrSeparatedValues fin(personFile);
    while (fin.nextline()) {
      ll person = fin.nextuint();
      fin.next();
      fin.next();
      fin.next();
      int y = fin.nextuint(), m = fin.nextuint(), d = fin.nextuint();
      tmpBirthday.push_back(make_pair(person, getBirthdayCode(y, m, d)));
      N++;
    }
    
    birthday.resize(N);
    for (int i = 0; i < N; i++){
      birthday[tmpBirthday[i].first] = tmpBirthday[i].second;
    }
  }
  
  {
    OrSeparatedValues fin(graphFile);
    while (fin.nextline()) {
      ll u = fin.nextuint();
      ll v = fin.nextuint();
      if (u < v){
        links.push_back(make_pair(u, v));
      }
    }
  }

  {
    OrSeparatedValues fin(tagListFile);
    while (fin.nextline()) {
      ll tag = fin.nextulong();
      tagToId[tag] = T;
      tagName.push_back(fin.next());
      T++;
    }
  }
    
  vertexTag     = vector<vector<LL> >(N);
  vertexTagSize = vector<int>(N + 1);
  
  {
    OrSeparatedValues fin(vertexTagFile);

    tagCount.resize(T);
    vector<pair<int, int> > tmpTagCount(T);
    vector<int> tagRank(T);
    
    for (int i = 0; i < T; i++){
      tmpTagCount[i] = make_pair(0, i);
    }
    
    while (fin.nextline()) {
      ll person  = fin.nextulong(), tag = fin.nextulong();
      int tagId = tagToId[tag];
      vertexTag[person].push_back(tagId);
      tmpTagCount[tagId].first--;
    }
    
    sort(tmpTagCount.begin(), tmpTagCount.end());
    
    for (int i = 0; i < T; i++){
      tagCount[i]                    = -tmpTagCount[i].first;
      tagRank[tmpTagCount[i].second] = i;
    }
    {
      vector<string> tmpTagName(T);
      for (int i = 0; i < T; i++){
        tmpTagName[tagRank[i]] = tagName[i];
      }
      tagName = tmpTagName;
    }
    
    for (int i = 0; i < N; i++){
      for (size_t j = 0; j < vertexTag[i].size(); j++){
        vertexTag[i][j] = tagRank[vertexTag[i][j]];
      }
    }
    
    rep(i, N){
      sort(vertexTag[i].begin(), vertexTag[i].end());
    }
    
    vertexTagSize[0] = 0;
    rep(i, N){
      vertexTagSize[i + 1] = vertexTagSize[i] + vertexTag[i].size();
    }
  }
  
  string line;
  ifstream ifsQuery(queryFile.c_str());
  while (getline(ifsQuery, line)) {
    if (line[5] != '2') continue;
    int k, y, m, d;
    sscanf(line.c_str(), "query2(%d,%d-%d-%d", &k, &y, &m, &d);
    queries.push_back(make_pair(getBirthdayCode(y, m, d), k));
    Q++;
  }
}
