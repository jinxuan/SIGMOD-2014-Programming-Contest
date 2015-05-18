#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <tuple>
#include "common.hpp"
#include "osv.hpp"
using namespace std;

const bool debug3 = false;

namespace {
  map<string, vector<int> > PlaceName2PlaceID;
  map<int, int> PlaceID2id;
  vector<vector<int> > PeopleIn;
  map<int, int> OrganizationID2PlaceID;
  vector<vector<int> > PersonGraph;
  vector<vector<int> > PlaceGraph;
  vector<vector<int> > Interests;
}

void collect_persons_dfs(int place, vector<int>& people){
  for (size_t i = 0; i < PeopleIn[place].size(); i++){
    people.push_back(PeopleIn[place][i]);
  }
  
  for (size_t i = 0; i < PlaceGraph[place].size(); i++){
    collect_persons_dfs(PlaceGraph[place][i], people);
  }
}

vector<int> collect_people_in(const vector<int> &places){
  vector<int> people;
  for (size_t i = 0; i < places.size(); i++){
    collect_persons_dfs(places[i], people);
  }
  sort(people.begin(), people.end());
  people.erase(unique(people.begin(), people.end()), people.end());
  return people;
}

int similarity(const vector<int> &interests_a, const vector<int> &interests_b){
  int    res   = 0;
  size_t pos_a = 0;
  size_t pos_b = 0;
  while(pos_a < interests_a.size() && pos_b < interests_b.size()){
    int i_a = interests_a[pos_a];
    int i_b = interests_b[pos_b];
    if(i_a == i_b){
      pos_a++;
      pos_b++;
      res++;
    }else if(i_a < i_b){
      pos_a++;
    } else {
      pos_b++;
    }
  }
  return res;
}

int similarity(const vector<bool> &interests_a, const vector<int> &interests_b){
  int res = 0;
  for (size_t i = 0; i < interests_b.size(); i++){
    res += interests_a[interests_b[i]];
  }
  return res;
}

bool reachable(int person1,
               int person2,
               vector<bool> &visit1,
               vector<bool> &visit2,
               int threshold){
  assert(person1 != person2);
  
  bool res = false;
  int dist = 0;
  int cur  = 0;
  
  queue<int> que1[2], que2[2];
  vector<int> updated1, updated2;
  
  que1[cur].push(person1);
  que2[cur].push(person2);
  updated1.push_back(person1);
  updated2.push_back(person2);
  visit1[person1] = true;
  visit2[person2] = true;
  
  while ((!que1[cur].empty() || !que2[cur].empty()) && dist < threshold){
    dist++;
    
    while (!que1[cur].empty() && dist <= threshold){
      int v = que1[cur].front(); que1[cur].pop();
      
      for (size_t i = 0; i < PersonGraph[v].size(); i++){
        int w = PersonGraph[v][i];
        
        if (visit2[w]){
          res = true;
          goto END;
        }
        
        if (!visit1[w]){
          visit1[w] = true;
          que1[1 - cur].push(w);
          updated1.push_back(w);
        }
      }
    }
    
    dist++;
    
    while (!que2[cur].empty() && dist <= threshold){
      int v = que2[cur].front(); que2[cur].pop();
      
      for (size_t i = 0; i < PersonGraph[v].size(); i++){
        
        int w = PersonGraph[v][i];
        
        if (visit1[w]){
          res = true;
          goto END;
        }
        
        if (!visit2[w]){
          visit2[w] = true;
          que2[1 - cur].push(w);
          updated2.push_back(w);
        }
      }
    }
    cur = 1 - cur;
  }
 END:
  for (size_t i = 0; i < updated1.size(); i++) visit1[updated1[i]] = false;
  for (size_t i = 0; i < updated2.size(); i++) visit2[updated2[i]] = false;
  return res;
}

void remove_redundant_people(vector<int> &people){
  map<int, int> interest_sum;
  for (size_t i = 0; i < people.size(); i++){
    vector<int> &I = Interests[people[i]];
    for (size_t j = 0; j < I.size(); j++){
      interest_sum[I[j]]++;
    }
  }
  
  int num_of_person = 0;
  for (size_t i = 0; i < people.size(); i++){
    bool flag = false;
    vector<int> &I = Interests[people[i]];
    
    for (size_t j = 0; j < I.size(); j++){
      if (interest_sum[I[j]] > 1) flag = true;
    }
    if (flag){
      people[num_of_person++] = people[i];
    }
  }

  people.resize(num_of_person);
}

int compress_interests(vector<vector<int> > &interests){
  map<int, int> interest_id;
  map<int, int> interest_count;

  for (size_t i = 0; i < interests.size(); i++){
    for (size_t j = 0; j < interests[i].size(); j++){
      interest_count[interests[i][j]]++;
    }
  }
  
  for (size_t i = 0; i < interests.size(); i++){
    int num_of_interests = 0;
    
    for (size_t j = 0; j < interests[i].size(); j++){
      if (interest_count[interests[i][j]] <= 1) continue;
      
      if (interest_id.find(interests[i][j]) == interest_id.end()){
        int new_id = interest_id.size();
        interest_id[interests[i][j]] = new_id;
      }
      interests[i][num_of_interests++] = interest_id[interests[i][j]];
    }
    interests[i].resize(num_of_interests);
    sort(interests[i].begin(), interests[i].end());
  }
  
  return interest_id.size();
}

vector<pair<int, int> > compute_lexicographic_pairs(const vector<int> &people,
                                                    size_t k,
                                                    int hop){
  vector<pair<int, int> > res;
  
  for (size_t i = 0; i < people.size(); i++){
    vector<int> dist(PersonGraph.size(), -1);
    dist[people[i]] = 0;
    queue<int> que;
    que.push(people[i]);
    while (!que.empty()){
      int v = que.front(); que.pop();
      for (auto it = PersonGraph[v].begin(); it != PersonGraph[v].end(); it++){
        if (dist[*it] == -1 && dist[v] + 1 <= hop){
          dist[*it] = dist[v] + 1;
          que.push(*it);
        }
      }
    }
    for (size_t j = i + 1; j < people.size(); j++){
      if (dist[people[j]] >= 0){
        res.push_back(make_pair(people[i], people[j]));
      }
    }
    if (res.size() >= k) break;
  }
  if (res.size() > k) res.resize(k);
  return res;
}

vector<pair<int,int> > query3(size_t k, int hop, const string& place){
  
  vector<int>  people(collect_people_in(PlaceName2PlaceID[place]));
  vector<bool> visit1(PersonGraph.size(), false);
  vector<bool> visit2(PersonGraph.size(), false);
  vector<pair<int, int> > topK;
  vector<pair<int, int> > lexicographic_pairs;
  lexicographic_pairs = compute_lexicographic_pairs(people, k, hop);
  
  remove_redundant_people(people);
  vector<vector<int> > interests;
  for (size_t i = 0; i < people.size(); i++)
    interests.push_back(Interests[people[i]]);
  int num_of_interests = compress_interests(interests);
  vector<vector<int> > person_has_interests(num_of_interests);

  for (size_t i = 0; i < people.size(); i++){
    for (size_t j = 0; j < interests[i].size(); j++){
      person_has_interests[interests[i][j]].push_back(i);
    }
  }


  
  priority_queue<tuple<int, int, int> > que;
  vector<set<int> > used(people.size());
  vector<bool> has_common_interests(people.size(), false);
  vector<bool> interests_of_person(num_of_interests, false);
  
  for (size_t i = 0; i < people.size(); i++){
    que.push(make_tuple<int, int, int>(interests[i].size(), -i, -i));
  }

  while (topK.size() < k && !que.empty()){
    auto pair = que.top(); que.pop();
    int id1     = -get<1>(pair);
    int id2     = -get<2>(pair);
    int person1 = people[id1];
    int person2 = people[id2];
    used[id1].insert(id2);
        
    if (person1 != person2 && reachable(person1, person2, visit1, visit2, hop)){
      topK.push_back(make_pair(person1, person2));
    }

    {
      int next = -1;
      int best = 0;
      
      vector<int> next_cands;
      for (size_t i = 0; i < interests[id1].size(); i++){
        const vector<int> &ps = person_has_interests[interests[id1][i]];
        
        size_t j = upper_bound(ps.begin(), ps.end(), id1) - ps.begin();
        for (; j < ps.size(); j++){
          if (!used[id1].count(ps[j]) && !has_common_interests[ps[j]]){
            next_cands.push_back(ps[j]);
            has_common_interests[ps[j]] = true;
          }
        }
        interests_of_person[interests[id1][i]] = true;
      }
      
      for (size_t i = 0; i < next_cands.size(); i++){
        int sim = similarity(interests_of_person, interests[next_cands[i]]);
        if (sim > best || (sim == best && next_cands[i] < next)){
          best = sim;
          next = next_cands[i];
        }
        has_common_interests[next_cands[i]] = false;
      }

      for (size_t i = 0; i < interests[id1].size(); i++){
        interests_of_person[interests[id1][i]] = false;
      }
      
      if (next != -1) que.push(make_tuple(best, -id1, -next));
    }
  }
  
  if (topK.size() < k){
    set<pair<int ,int> > topK_pair;
    for (size_t i = 0; i < topK.size(); i++) topK_pair.insert(topK[i]);
    for (size_t i = 0; i < lexicographic_pairs.size() && topK.size() < k; i++){
      if (topK_pair.count(lexicographic_pairs[i])) continue;
      topK.push_back(lexicographic_pairs[i]);
    }
  }

  return topK;
}

void load_query(const string &query_file, vector<tuple<int, int, string> > &qs){
  ifstream ifs(query_file.c_str());
  string   in;
  
  while(getline(ifs, in)){
    if(in[5] != '3') continue;
    for (size_t i = 0; i < in.size(); i++){
      if(in[i] == ',') in[i] = ' ';
      if(in[i] == '(') in[i] = ' ';
      if(in[i] == ')') in[i] = ' ';
    }
      
    istringstream ss(in);
    int k, h;
    string p, queryname;
    ss >> queryname >> k >> h >> p;
    qs.push_back(make_tuple(k, h, p));
  }
  ifs.close();
}

int main(int argc, char *argv[])
{
  double start = get_current_time_sec();
  
  if (argc != 3){
    cerr << "Usage: " << argv[0]
         << " (path_to_data_directory)"
         << " (path_to_query_file)" << endl;
    exit(EXIT_FAILURE);
  }
  
  int V = 0;
  const string data_dir   = argv[1] + string("/");
  const string query_file = argv[2];
  
  {
    double s = get_current_time_sec();
    vector<pair<LL, LL> > es;
    load_pairs(data_dir + "person_knows_person.csv", es);
    for (auto it = es.begin(); it != es.end(); it++){
      V = max(V, (int)max(it->first, it->second) + 1);
    }
    
    PersonGraph.resize(V);
    for (auto it = es.begin(); it != es.end(); it++){
      PersonGraph[it->first].push_back(it->second);
    }
    message(3, "Load person-graph", s, debug3);
  }
  
  
  {
    double s = get_current_time_sec();
    OrSeparatedValues osv(data_dir + "place.csv");
    while(osv.nextline()){
      int    placeID   = osv.nextuint();
      string placeName = osv.next();
      int cnt = PlaceID2id.size();
      PlaceID2id[placeID] = cnt;
      PlaceName2PlaceID[placeName].push_back(cnt);
    }
    message(3, "Load places", s, debug3);
  }
  
  PlaceGraph.resize(PlaceID2id.size());
  
  {
    double s = get_current_time_sec();
    OrSeparatedValues osv(data_dir + "place_isPartOf_place.csv");
    while(osv.nextline()){
      int placeA = osv.nextuint();
      int placeB = osv.nextuint();
      PlaceGraph[PlaceID2id[placeB]].push_back(PlaceID2id[placeA]);
    }
    message(3, "Load place-graph", s, debug3);
  }
  
  {
    double s = get_current_time_sec();
    OrSeparatedValues osv(data_dir + "organisation_isLocatedIn_place.csv");
    while(osv.nextline()){
        int organizationID = osv.nextuint();
        int placeID        = osv.nextuint();
        OrganizationID2PlaceID[organizationID] = placeID;      
    }
    message(3, "Load org-place", s, debug3);
  }
  
  PeopleIn.resize(PlaceID2id.size());
  {
    double s = get_current_time_sec();
    {
      OrSeparatedValues osv(data_dir + "person_isLocatedIn_place.csv");
      while(osv.nextline()){
        int personID = osv.nextuint();
        int placeID  = osv.nextuint();
        int placeid  = PlaceID2id[placeID];
        PeopleIn[placeid].push_back(personID);
      }
    }
    
    {
      OrSeparatedValues osv(data_dir + "person_workAt_organisation.csv");
      while(osv.nextline()){
        int personID = osv.nextuint();
        int organizationID = osv.nextuint();
        int PlaceID = OrganizationID2PlaceID[organizationID];
        int placeid = PlaceID2id[PlaceID];      
        PeopleIn[placeid].push_back(personID);
      }
    }
  
    {
      OrSeparatedValues osv(data_dir + "person_studyAt_organisation.csv");
      while(osv.nextline()){
        int personID       = osv.nextuint();
        int organizationID = osv.nextuint();
        int PlaceID = OrganizationID2PlaceID[organizationID];
        PeopleIn[PlaceID2id[PlaceID]].push_back(personID);
      }
    }
    message(3, "Load place-people", s, debug3);
  }
  
  Interests.resize(V);
  
  {
    double s = get_current_time_sec();
    OrSeparatedValues osv(data_dir + "person_hasInterest_tag.csv");
    while(osv.nextline()){
      int personID   = osv.nextuint();
      int interestID = osv.nextuint();
      Interests[personID].push_back(interestID);
    }

    
    message(3, "Load interests", s, debug3);
  }

  
  vector<tuple<int, int, string> > query_input;
  {
    double   s = get_current_time_sec();
    load_query(query_file, query_input);
    message(3, "Load query", s, debug3);
  }
  
  vector<vector<pair<int, int> > > query_output(query_input.size());
  
  {
    double s = get_current_time_sec();
    for (size_t i = 0; i < query_input.size(); i++){
      auto query      = query_input[i];
      query_output[i] = query3(get<0>(query), get<1>(query), get<2>(query));
    }
    message(3, "Queries", s, debug3);
  }
  
  for(size_t i = 0; i < query_output.size(); i++){
    auto &ans = query_output[i];
    for (size_t i = 0; i < ans.size(); i++){
      cout << ans[i].first << "|" << ans[i].second
           << (i + 1 == ans.size() ? "\n" : " " );
    }
  }
  
  message(3, "Total", start, true);
  return 0;
}
