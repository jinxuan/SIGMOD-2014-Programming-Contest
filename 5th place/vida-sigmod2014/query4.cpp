//----------------------------------------------------------------------------
//
// Copyright (C) 2014, NYU Polytechnic School of Engineering. 
// All rights reserved.
// Contact: fchirigati@nyu.edu, kien.pham@nyu.edu,
//          tuananh@nyu.edu, huy.vo@nyu.edu
//
// This file is part of the VIDA team submission for the SIGMOD 2014
// Programming Contest.
//
// "Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.
//  - Redistributions in binary form must reproduce the above copyright 
//    notice, this list of conditions and the following disclaimer in the 
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of NYU nor the names of its 
//    contributors may be used to endorse or promote products derived from 
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
//
//----------------------------------------------------------------------------

#include "query.hpp"
#include <deque>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

struct Forum
{
  typedef std::vector<uint64_t> IdVector;
  typedef boost::unordered_map<uint64_t, uint32_t> IdMap;

  IdMap    idMap;

  Forum(const char *dataPath=0)
  {
    if (dataPath) {
      this->readForum(dataPath);
    }
  }

  void readForum(const char *dataPath)
  {
    IdVector id;
    id.reserve(SharedData::instance()->person->size()*16);
    {
      std::string fn = std::string(dataPath) + "/forum.csv";
      boost::iostreams::mapped_file_source file(fn);
      const char *cur = file.data();
      const char *end = cur + file.size();
      while (*cur!='\n') cur++;
      while (++cur<end) {
        uint64_t p = 0;
        while (*cur>='0' && *cur<='9')
          p = p*10 + (*cur++-'0');
        id.push_back(p);
        while (*cur!='\n') cur++;
      }
    }
    this->idMap.clear();
    this->idMap.rehash(id.size()/this->idMap.max_load_factor()+1);
    for (uint32_t i=0; i<id.size(); i++)
      this->idMap[id[i]] = i;
  }

  inline size_t size() const
  {
    return this->idMap.size();
  }

  inline uint32_t normalized(uint64_t p) const
  {
    IdMap::const_iterator it = this->idMap.find(p);
    if (it!=this->idMap.end())
      return it->second;
    return -1;
  }

};

struct TagPerson
{
  typedef boost::unordered_map<uint64_t, uint32_t> ForumIdMap;
  typedef std::pair<uint64_t, uint32_t> ForumPair;
  typedef std::vector<ForumPair>        ForumPairVector;
  typedef std::vector<uint32_t>         IdVector;
  typedef std::pair<uint32_t,uint32_t>  IdPair;
  typedef std::vector<IdPair>           IdPairVector;

  Forum * forum;
  Tag * tag;
  const Person * person;
  IdVector tagForum, forumPerson;
  IdPairVector tagForumIndex;
  IdPairVector forumPersonIndex;

  TagPerson(Forum * forum, Tag *tag, const Person *person, const char *dataPath=0)
  {
    this->forum = forum;
    this->tag = tag;
    this->person = person;
    if (dataPath) {
      // boost::thread worker = boost::thread(&TagPerson::readTagForum, this, dataPath);
      this->readTagForum(dataPath);
      this->readForumPerson(dataPath);
      // worker.join();
    }
  }

  void readTagForum(const char * dataPath)
  {
    {
      IdPairVector tag_to_forum;
      std::string fn = std::string(dataPath) + "/forum_hasTag_tag.csv";
      boost::iostreams::mapped_file_source file(fn);
      const char *cur = file.data();
      const char *end = cur + file.size();
      uint64_t p = 0, fid = 0;
      while (*cur!='\n') cur++;
      while (++cur<end) {
        if (*cur>='0' && *cur<='9')
          p = p*10 + (*cur-'0');
        else {
          if (*cur=='|') {
            fid = p;
          }
          else if (*cur=='\n')
            tag_to_forum.push_back(IdPair(this->tag->normalized(p),
                                          this->forum->normalized(fid)));
          else
            exit(2);
          p = 0;
        }
      }
      hsort::sort(tag_to_forum.begin(), tag_to_forum.end()); 
      this->tagForum.resize(tag_to_forum.size());
      this->tagForumIndex.resize(this->tag->size(), IdPair(0, 0));
      IdPair last(0, 0);
      for (unsigned i=0; i<tag_to_forum.size(); i++) {
        uint32_t tid = tag_to_forum[i].first;
        if (tid!=last.first) {
          this->tagForumIndex[last.first] = IdPair(last.second, i-last.second);
          last = IdPair(tid, i);
        }
        this->tagForum[i] = tag_to_forum[i].second;
      }
      this->tagForumIndex[last.first] = IdPair(last.second, tag_to_forum.size()-last.second);
    }
  }

  void readForumPerson(const char * dataPath)
  {
    this->forumPerson.clear();
    this->forumPersonIndex.resize(this->forum->size());
    {
      std::string fn = std::string(dataPath) + "/forum_hasMember_person.csv";
      boost::iostreams::mapped_file_source file(fn);
      const char *cur = file.data();
      const char *end = cur + file.size();
      while (*cur!='\n') cur++;
      this->forumPerson.reserve(Comment::estimateLineCount(cur, end, file.size()));
      uint64_t last_fid=-1;
      IdPair lastIndex(0, 0);
      while (++cur<end) {
        uint64_t p=0;
        for (; *cur>='0' && *cur<='9'; ++cur)
          p = p*10 + (*cur-'0');
        if (p!=last_fid) {
          this->forumPersonIndex[lastIndex.first] = IdPair(lastIndex.second, this->forumPerson.size()-lastIndex.second);
          last_fid = p;
          lastIndex = IdPair(this->forum->normalized(last_fid), this->forumPerson.size());
        }

        p = 0;
        for (cur++; *cur>='0' && *cur<='9'; ++cur)
          p = p*10 + (*cur-'0');
        this->forumPerson.push_back(this->person->normalized(p));
        
        while (*cur!='\n') cur++;
      }
      this->forumPersonIndex[lastIndex.first] = IdPair(lastIndex.second, this->forumPerson.size()-lastIndex.second);
    }
  }

  uint32_t getPersonByTag(const char * t, BitSet &personMap) const
  {
    uint32_t count = 0;
    uint32_t tid = this->tag->tagNameMap[t];
    uint32_t min_tag = this->tagForumIndex[tid].first;
    uint32_t max_tag = min_tag+this->tagForumIndex[tid].second;
    for (uint32_t i=min_tag; i<max_tag; i++) {
      uint32_t min_forum = this->forumPersonIndex[this->tagForum[i]].first;
      uint32_t max_forum = min_forum + this->forumPersonIndex[this->tagForum[i]].second;
      for (uint32_t i_forum=min_forum; i_forum<max_forum; i_forum++) {
        if (personMap[this->forumPerson[i_forum]]==0) {
          personMap[this->forumPerson[i_forum]] = 1;
          count++;
        }
      }
    }
    return count;
  }
};

struct Query4
{
  typedef std::pair<uint32_t, uint32_t> IdPair;
  typedef std::vector<uint32_t>         U32Vector;
  typedef std::vector<IdPair> IdPairVector;
  typedef boost::unordered_map<IdPair, uint32_t> IdPairMap;
  typedef boost::unordered_map<uint64_t, uint32_t> U64Map;
  typedef std::pair<double, uint64_t> TopK;
  typedef std::vector<TopK> TopKVector;

  const TagPerson *tp;
  PersonGraph     *pg;
  uint32_t         personSize;
  U32Vector        QR;
  uint32_t        *R;
  uint32_t        *D;
  uint32_t        *E;
  
  boost::mutex     pMutex;
  uint32_t         pCurrent;
  uint32_t         pMax;

  Query4(PersonGraph *G, const TagPerson *tp)
  {
    this->pg = G;
    this->tp = tp;
    this->personSize = this->tp->person->size();
    this->QR.resize(3*this->personSize+1 + this->pg->adjacencyList.size() + 8*this->personSize);
    this->R = &this->QR[0];
    this->D = this->R + personSize;
    this->E = this->D + personSize;
  }

  struct topk_comparator {
    bool operator()(const TopK &lhs, const TopK &rhs) const {
      if (lhs.first == rhs.first)
        return (lhs.second > rhs.second);
      return (lhs.first < rhs.first);
    }
  };
  typedef std::set<TopK, Query4::topk_comparator> TopKSet;

  template<typename IdType>
  void createAdjacencyList(uint32_t vertexCount, uint32_t *vMap, BitSet &personMap, uint32_t *E, uint32_t *A)
  {
    uint32_t edgeCount = 0;
    IdType *edges = (IdType*)(A);
    for (int nu=0; nu<vertexCount; nu++) {
      uint32_t u = this->D[nu];
      E[nu] = edgeCount;
      for (int k=this->pg->edgeIndex[u]; k<this->pg->edgeIndex[u+1]; k++) {
        uint32_t v = this->pg->adjacencyList[k];
        if (personMap[v]) {
          edges[edgeCount++] = vMap[v];
        }
      }
    }
    E[vertexCount] = edgeCount;
  }

  std::string compute(int k, const char * t)
  {
    BitSet personMap(this->personSize);
    uint32_t n = this->tp->getPersonByTag(t, personMap);
    if (n==0) return "";
    this->D     = this->R + n;
    this->E     = this->D + n;
    uint32_t *A = this->E + n+1;
    uint32_t *Q = A + this->pg->adjacencyList.size();
    
    bool use16bit = n<65536;
    {
      for (int u=0,cnt=0; u<this->personSize; u++) {
        if (personMap[u]) {
          this->D[cnt] = u;
          Q[u] = cnt++;
        }
      }
    }
    if (use16bit)
      this->createAdjacencyList<uint16_t>(n, Q, personMap, E, A);
    else 
      this->createAdjacencyList<uint32_t>(n, Q, personMap, E, A);
    Q = A + E[n];
 
    IdPairVector rank(n);
    {
      personMap.resize(n);
      personMap.reset();
      for (int i=0; i<n; i++) {
        if (!personMap[i]) {
          if (use16bit)
            this->compute_r_p(i, personMap, n, E, (uint16_t*)A, Q);
          else
            this->compute_r_p(i, personMap, n, E, A, Q);
        }
        uint32_t deg = E[i+1]-E[i];
        rank[i] = IdPair(~deg, i);
      }
      hsort::sort(rank.begin(), rank.end());
    }
    k = std::min(k, (int)n);
    uint32_t topPos = k + (n-k) % NBITS;
    double top_min = -1;
    const int nthread = 8;
    TopKSet pc[nthread];
    boost::thread_group workers;
    for (int i=0; i<nthread; i++)
      if (use16bit)
        workers.create_thread(boost::bind(&Query4::call<uint16_t>, this, IdPair(i, topPos), nthread, top_min,
                                          IdPair(n, k), boost::cref(rank), (uint16_t*)A, Q + i*n, boost::ref(pc[i])));
      else
        workers.create_thread(boost::bind(&Query4::call<uint32_t>, this, IdPair(i, topPos), nthread, top_min,
                                          IdPair(n, k), boost::cref(rank), A, Q + i*n, boost::ref(pc[i])));
    workers.join_all();
 
    for (int i=0; i<nthread; i++)
      if (pc[i].size()>0 && (top_min<0 || top_min>pc[i].begin()->first))
        top_min = pc[i].begin()->first;

    this->pCurrent = topPos;
    this->pMax = n;
    for (int i=0; i<nthread; i++)
      if (use16bit)
        workers.create_thread(boost::bind(&Query4::call64<uint16_t>, this, nthread, top_min,
                                          IdPair(n, k), boost::cref(rank), (uint16_t*)A, boost::ref(pc[i])));
      else
        workers.create_thread(boost::bind(&Query4::call64<uint32_t>, this, nthread, top_min,
                                          IdPair(n, k), boost::cref(rank), A, boost::ref(pc[i])));
    workers.join_all();

    TopKVector personCentrality;
    for (int i=0; i<nthread; i++)
      personCentrality.insert(personCentrality.end(), pc[i].begin(), pc[i].end());

    std::sort(personCentrality.begin(), personCentrality.end(), topk_comparator());
      
    std::string results = "";
    char result[64];
    for (TopKVector::reverse_iterator it=personCentrality.rbegin(); k>0 && it!=personCentrality.rend(); it++, --k) {
      sprintf(result, "%llu ", it->second);
      results += result;
    }
    if (!results.empty())
      results.resize(results.size()-1);
    return results;
  }
  
#define EXPAND(u) {                                     \
    if (!curF[u]) continue;                             \
    register BitMask bm = curF[u] & active;             \
    curF[u] = 0;                                        \
    if (!bm) continue;                                  \
    for(int k=this->E[u],e=this->E[u+1]; k<e; k++) {    \
      register uint32_t v = adjacencyList[k];           \
      register BitMask nextV = bm & (~VV[v]);           \
      if (nextV) {                                      \
        nextF[v] |= nextV;                              \
        VV[v] |= nextV;                                 \
        do {                                            \
          --r_p[__builtin_ctzl(nextV)];                 \
        } while (nextV &= (nextV-1));                   \
      }                                                 \
    }                                                   \
  }

#define UPDATE() {                                      \
    for (BitMask a=active; a; a&=(a-1)) {               \
      register int j = __builtin_ctzl(a);               \
      register BitMask nshift = ~(1UL << j);            \
      if (r_p[j]==0)                                    \
        active &= nshift;                               \
      else {                                            \
        s_p[j] += r_p[j];                               \
        if (s_p[j]>s_p_max[j]) {                        \
          s_p[j] = 0;                                   \
          active &= nshift;                             \
        }                                               \
      }                                                 \
    }                                                   \
  }

  template<typename IdType>
  void call64(int stride, double current_min, IdPair nk, const IdPairVector &rank, IdType *adjacencyList, TopKSet &personCentrality)
  {
    BitMask *F = (BitMask*)malloc(3*nk.first*sizeof(BitMask));
    BitMask *FF[2] = {F, F + nk.first};
    BitMask *VV    = F + 2*nk.first;

    uint32_t i=0;
    while (1) {
      {
        this->pMutex.lock();
        i = this->pCurrent;
        this->pCurrent += NBITS;
        this->pMutex.unlock();
      }
      if (i>=this->pMax) break;

      // initialize
      memset(F, 0, 3*nk.first*sizeof(BitMask));
      IdType   r_p[NBITS];
      uint32_t s_p[NBITS];
      uint32_t s_p_max[NBITS];
      double value[NBITS];
      BitMask *curF=FF[0], *nextF=FF[0];
      BitMask active = -1;
      for (BitMask j=0,shift=1; j<NBITS; j++, shift<<=1) {
        uint32_t p = rank.at(i+j).second;
        r_p[j] = this->R[p]-1;
        value[j] = (nk.first<2)?0:(double)r_p[j]*r_p[j]/(double)(nk.first-1);
        s_p[j] = r_p[j];
        s_p_max[j] = current_min==0?INT_MAX:value[j]/current_min;
        if (r_p[j]<2) {
          active &= ~shift;
          s_p[j] = 0;
        }
        VV[p] |= shift;
        nextF[p] |= shift;
      }

      // performing bfs for hop=1
      for (int hop=1; hop<2 && active; hop++, curF=nextF) {
        nextF = FF[hop%2];
        for (uint32_t j=0; j<NBITS; j++) {
          uint32_t u = rank.at(i+j).second;
          EXPAND(u);
        }
        UPDATE();
      }

      // performing bfs for the rest
      for (int hop=2; active; hop++, curF=nextF) {
        nextF = FF[hop%2];
        for (uint32_t u=0; u<nk.first; u++)
          EXPAND(u);
        UPDATE();
      }

      // finalize
      for (uint32_t j=0; j<NBITS; j++) {
        uint32_t p = rank.at(i+j).second;
        double cc = s_p[j]?(value[j]/s_p[j]):0;
        if (cc>=current_min) {
          personCentrality.insert(std::make_pair(cc, this->tp->person->denormalized(this->D[p])));
          if (personCentrality.size()>nk.second) {
            personCentrality.erase(personCentrality.begin());
            current_min = personCentrality.begin()->first;
          }
        }
      }
    }
    free(F);
  }

  template<typename IdType>
  void compute_r_p(uint32_t pid, BitSet &personMap, uint32_t n, uint32_t *edgeIndex, IdType *adjacencyList, uint32_t *Q)
  {
    int front=0, count=0;
    personMap[pid] = 1;
    Q[count++] = pid;
    while (front<count) {
      uint32_t u = Q[front++];
      for(int k=edgeIndex[u]; k<edgeIndex[u+1]; k++) {
        uint32_t id = adjacencyList[k];
        if (!personMap[id]) {
          personMap[id] = 1;
          Q[count++] = id;
        }
      }
    }
    for (int i=0; i<count; i++)
      this->R[Q[i]] = count;
  }  

  template<typename IdType>
  double closeness_centrality(uint32_t pid, double kmin, uint32_t n, uint32_t *edgeIndex, IdType *adjacencyList, uint32_t *Q)
  {
    uint32_t r_p = this->R[pid];
    if (r_p==1) return 0.0;
    uint32_t s_p = 0;
    double value = (n<1)?0:(double)(r_p-1)*(r_p-1)/(double)(n-1);
    uint32_t s_p_max = kmin==0?INT_MAX:value/kmin;    
    int front=0, count=0;
    BitSet V(n);
    V[pid] = 1;
    Q[count++] = pid;
    for (int hop=1; front<count; hop++) {
      if (s_p+(r_p-count)*hop>s_p_max)
        return 0.0;
      for (int end=count; front<end; ++front) {
        uint32_t u = Q[front];
        for(int k=edgeIndex[u]; k<edgeIndex[u+1]; k++) {
          uint32_t id = adjacencyList[k];
          if (!V[id]) {
            V[id] = 1;
            Q[count++] = id;
          }
        }
      }
      s_p += hop*(count-front);
    }
    if (s_p==0)
      return 0.0;
    return value/s_p;
  }

  template<typename IdType>
  void call(IdPair range, int stride, double current_min, IdPair nk, const IdPairVector &rank, IdType *adjacencyList, uint32_t *Q, TopKSet &personCentrality)
  {
    for (uint32_t i=range.first; i<range.second; i+=stride) {
      uint32_t p = rank.at(i).second;
      double cc = this->closeness_centrality(p, current_min, nk.first, this->E, adjacencyList, Q);
      if (cc>=current_min) {
        personCentrality.insert(std::make_pair(cc, this->tp->person->denormalized(this->D[p])));
        if (personCentrality.size()>nk.second) {
          personCentrality.erase(personCentrality.begin());
          current_min = personCentrality.begin()->first;
        }
      }
    }
  }

  static void solve(PersonGraph *G, const TagPerson *tp, std::deque<std::string*> &q4)
  {
    static boost::mutex mutex;
    Query4 solver(G, tp);
    while (true) {
      std::string *nextQ = 0;
      {
        boost::mutex::scoped_lock lock(mutex);
        if (!q4.empty()) {
          nextQ = q4.front();
          q4.pop_front();
        }
      }
      if (!nextQ) break;
      int k;
      char t[128] = {};
      sscanf(nextQ->c_str(), "query4(%d, %[^)])\n", &k, t);
      *nextQ = solver.compute(k, t);
    }
  }
};

struct Query4Data {
  Forum     * forum;
  TagPerson * tp;
  Query4    * solver;
};

void Query4Thread::initialize()
{
  if (this->queries->empty()) return;
  Query4Data *data = new Query4Data();
  data->forum  = new Forum(this->dataPath);
  data->tp     = new TagPerson(data->forum, SharedData::instance()->tag, SharedData::instance()->person, this->dataPath);
  data->solver = new Query4(SharedData::instance()->personGraph, data->tp);
  this->data = data;
}

void Query4Thread::perform()
{
  if (this->queries->empty()) return;
  this->initThread->join();
  Query4Data *data = (Query4Data*)this->data;
  for (unsigned i=0; i<this->queries->size(); i++) {
    int k;
    char t[128] = {};
    sscanf(this->queries->at(i)->c_str(), "query4(%d, %[^)])\n", &k, t);
    this->queries->at(i)->assign(data->solver->compute(k, t));
  }
  delete data->solver;
  delete data->tp;
  delete data->forum;
  delete data;
}
