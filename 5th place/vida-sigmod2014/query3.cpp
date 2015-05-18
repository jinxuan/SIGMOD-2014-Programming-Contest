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

struct Place
{
  typedef std::vector<uint64_t>                                      IdVector;
  typedef boost::unordered_map<uint64_t, uint32_t>                   IdMap;
  typedef boost::unordered_map<std::string, uint32_t>                PlaceMap;
  typedef boost::unordered_map< std::string, std::vector<uint32_t> > PlaceVectorMap;
  typedef std::pair<uint64_t, std::string>                           IdPair;
  typedef std::vector<IdPair>                                        IdPairVector;

  IdVector       id;
  IdMap          idMap;
  bool           defaultId;
  PlaceMap       placeMap;
  PlaceVectorMap placeVectorMap;

  Place(const char *dataPath=0)
  {
    if (dataPath) {
      this->readPlace(dataPath);
    }
  }

  void readPlace(const char *dataPath)
  {
    IdPairVector pid_to_name;

    this->id.clear();
    this->idMap.clear();
    std::string fn = std::string(dataPath) + "/place.csv";
    FILE *fi = fopen(fn.c_str(), "r");
    fscanf(fi, "%*[^\n]\n");
    while (!feof(fi)) {
      uint64_t pid;
      char char_name[128] = {};
      fscanf(fi, "%llu|%[^|]|%*[^\n]\n", &pid, char_name);
      pid_to_name.push_back(IdPair(pid, std::string(char_name)));
    }
    fclose(fi);
    std::sort(pid_to_name.begin(), pid_to_name.end());
    for (uint32_t i=0; i<pid_to_name.size(); i++) {
      this->id.push_back(pid_to_name[i].first);
      this->idMap[this->id[i]] = i;

      std::string name = pid_to_name[i].second;
      PlaceMap::iterator it = this->placeMap.find(name);
      if (it!=this->placeMap.end()) {
        if (it->second==0)
          this->placeVectorMap[name].push_back(i);
        else {
          std::vector<uint32_t> & place_id = this->placeVectorMap[name];
          place_id.push_back(it->second - 1);
          place_id.push_back(i);
          it->second = 0;
        }
      } else {
        this->placeMap[name] = i + 1;
      }

    }
    this->defaultId = (this->id.back()+1)==this->id.size();
  }

  inline size_t size() const
  {
    return this->id.size();
  }

  inline uint32_t normalized(uint64_t p) const
  {
    if (this->defaultId)
      return p;
    IdMap::const_iterator it = this->idMap.find(p);
    if (it!=this->idMap.end())
      return it->second;
    return -1;
  }
};

struct PlaceGraph
{
  typedef std::pair<uint32_t, uint32_t> IdPair;
  typedef std::vector<IdPair>           IdPairVector;
  typedef std::vector<uint32_t>         IdVector;

  Place       * place;
  size_t        placeSize;

  IdVector             adjacencyList;
  std::vector<int>     edgeIndex;

  PlaceGraph(Place *P, const char *dataPath=0)
  {
    this->place = P;
    this->placeSize = this->place->size();
    if (dataPath) {
      this->buildGraph(dataPath);
    }
  }

  void buildGraph(const char *dataPath)
  {
    IdPairVector edges;

    std::string fn = std::string(dataPath) + "/place_isPartOf_place.csv";
    FILE *fi = fopen(fn.c_str(), "r");
    fscanf(fi, "%*[^\n]\n");
    while (!feof(fi)) {
      uint64_t p[2];
      fscanf(fi, "%llu|%llu\n", p, p+1);
      IdPair ip(this->place->normalized(p[1]), this->place->normalized(p[0]));
      edges.push_back(ip);
    }
    fclose(fi);

    hsort::sort(edges.begin(), edges.end());

    this->edgeIndex.clear();
    this->edgeIndex.resize(this->placeSize+1, 0);
    this->adjacencyList.resize(edges.size(), 0);

    for (unsigned i=0; i<edges.size(); i++)
      this->edgeIndex[edges[i].first+1]++;

    for (unsigned i=0; i<this->placeSize; i++)
      this->edgeIndex[i+1] += this->edgeIndex[i];

    for (unsigned i=0; i<edges.size(); i++)
      this->adjacencyList[i] = edges[i].second;

    this->V.resize(this->placeSize);
  }

  mutable BitSet V;
  void visitedPlaces(uint32_t src, IdVector &Q) const
  {
    if (this->V[src]) return;
    this->V[src] = 1;
    uint32_t front = Q.size();
    Q.push_back(src);
    while (front<Q.size()) {
      uint32_t u = Q[front++];
      for (int k=this->edgeIndex[u]; k<this->edgeIndex[u+1]; k++) {
        const uint32_t &id = this->adjacencyList[k];
        if (!this->V[id]) {
          this->V[id] = 1;
          Q.push_back(id);
        }
      }
    }
  }

  inline void getID(const char *name, IdVector &places) const
  {
    this->V.reset();
    std::string str_name(name);
    Place::PlaceMap::const_iterator it = this->place->placeMap.find(str_name);
    if (it==this->place->placeMap.end()) return;
    if (it->second==0) {
      const std::vector<uint32_t> &place_id = this->place->placeVectorMap.at(str_name);
      for (unsigned i=0; i<place_id.size(); i++)
        this->visitedPlaces(place_id[i], places);
    } else {
      this->visitedPlaces(it->second-1, places);
    }
  }  
};

struct PlaceHasPerson
{
  typedef std::pair<uint32_t, uint32_t>            IdPair;
  typedef std::vector<IdPair>                      IdPairVector;
  typedef boost::unordered_map<uint64_t, uint32_t> IdMap;
  typedef boost::unordered_set<uint32_t>           IdSet;
  typedef std::vector<IdSet>                       IdSetVector;

  const Place  *place;
  const Person *person;
  size_t        placeSize;
  size_t        personSize;

  IdSetVector placeHasPerson;

  PlaceHasPerson(const Person *person, const Place *place, const char* dataPath=0)
  {
    this->person = person;
    this->place = place;
    this->personSize = this->person->size();
    this->placeSize = this->place->size();
    if (dataPath) {
      this->buildVector(dataPath);
    }
  }

  void buildVector(const char* dataPath) {

    IdPairVector placeid_has_pid;
    IdMap oid_in_placeid;

    // Reading person_isLocatedIn_place
    std::string fn = std::string(dataPath) + "/person_isLocatedIn_place.csv";
    FILE *fi = fopen(fn.c_str(), "r");
    fscanf(fi, "%*[^\n]\n");
    while (!feof(fi)) {
      uint64_t person_id, place_id;
      fscanf(fi, "%llu|%llu\n", &person_id, &place_id);
      IdPair ip(this->place->normalized(place_id), this->person->normalized(person_id));
      placeid_has_pid.push_back(ip);
    }
    fclose(fi);

    // Reading organisation_isLocatedIn_place
    fn = std::string(dataPath) + "/organisation_isLocatedIn_place.csv";
    fi = fopen(fn.c_str(), "r");
    fscanf(fi, "%*[^\n]\n");
    while (!feof(fi)) {
      uint64_t oid, place_id;
      fscanf(fi, "%llu|%llu\n", &oid, &place_id);
      oid_in_placeid[oid] = this->place->normalized(place_id);
    }
    fclose(fi);

    // Reading person_workAt_organisation
    fn = std::string(dataPath) + "/person_workAt_organisation.csv";
    fi = fopen(fn.c_str(), "r");
    fscanf(fi, "%*[^\n]\n");
    while (!feof(fi)) {
      uint64_t person_id, oid;
      fscanf(fi, "%llu|%llu|%*[^\n]\n", &person_id, &oid);
      IdPair ip(oid_in_placeid[oid], this->person->normalized(person_id));
      placeid_has_pid.push_back(ip);
    }
    fclose(fi);

    // Reading person_studyAt_organisation
    fn = std::string(dataPath) + "/person_studyAt_organisation.csv";
    fi = fopen(fn.c_str(), "r");
    fscanf(fi, "%*[^\n]\n");
    while (!feof(fi)) {
      uint64_t person_id, oid;
      fscanf(fi, "%llu|%llu|%*[^\n]\n", &person_id, &oid);
      IdPair ip(oid_in_placeid[oid], this->person->normalized(person_id));
      placeid_has_pid.push_back(ip);
    }
    fclose(fi);

    hsort::sort(placeid_has_pid.begin(), placeid_has_pid.end());

    this->placeHasPerson.resize(this->placeSize);

    uint32_t current_placeid = placeid_has_pid.front().first;
    IdSet current_pid;
    current_pid.insert(placeid_has_pid.front().second);
    for(IdPairVector::iterator it = placeid_has_pid.begin(); it != placeid_has_pid.end(); it++)
      {
        if (it->first == current_placeid)
          current_pid.insert(it->second);
        else
          {
            this->placeHasPerson[current_placeid] = current_pid;
            current_placeid = it->first;
            current_pid.clear();
            current_pid.insert(it->second);
          }
      }
    this->placeHasPerson[current_placeid] = current_pid;

  }
};

struct TopK
{
  typedef std::vector<uint32_t>                   IdVector;
  typedef std::pair<uint32_t, uint32_t>           IdPair32;
  typedef std::vector<IdPair32>                   IdPairVector32;

  typedef std::pair<uint64_t, uint64_t>           IdPair;
  typedef std::pair<uint32_t, IdPair>             IdPairTopK;

  typedef std::vector<IdPairTopK> TopKVector;

  const PersonGraph    *personGraph;
  const PlaceGraph     *placeGraph;
  const PlaceHasPerson *placeHasPerson;
  Tag                  *personHasInterestTag;
  boost::mutex          pMutex;
  uint32_t              pCurrent;
  std::vector<uint32_t> masterQ;

  TopK(const PersonGraph    *personGraph,
       const PlaceGraph     *placeGraph,
       const PlaceHasPerson *placeHasPerson,
       Tag                  *personHasInterestTag)
  {
    this->personGraph          = personGraph;
    this->placeGraph           = placeGraph;
    this->placeHasPerson       = placeHasPerson;
    this->personHasInterestTag = personHasInterestTag;
  }

  std::string compute(int k, int h, const char *placeName)
  {
    IdVector places;
    places.reserve(this->placeGraph->placeSize);
    this->placeGraph->getID(placeName, places);

    BitSet persons(this->personGraph->personSize);
    int nPerson = 0;
    for (IdVector::iterator it = places.begin(); it!= places.end(); it++)
      {
        PlaceHasPerson::IdSet persons_set = this->placeHasPerson->placeHasPerson[*it];
        for (PlaceHasPerson::IdSet::iterator its = persons_set.begin(); its != persons_set.end(); its++) {
          if (!persons[*its]) {
            persons[*its] = 1;
            nPerson++;
          }
        }
      }


    if (nPerson==0) return "";

    //initialize priority queue
    const int nthread = 8;
    std::vector<TopKVector> topk_queue(nthread);
    uint32_t current_minimum = 0;
    IdPair32 *bp = &this->personHasInterestTag->byPerson[0];
    int bps = this->personHasInterestTag->byPerson.size();
    IdPairVector32 idp;
    idp.reserve(bps);
    for (int i=0; i<bps; i++)
      if (persons[bp[i].second])
        idp.push_back(bp[i]);
    bps = idp.size();

    bp = &idp[0];
    IdVector Q((~idp[0].first)+2, bps);
    for (int i=0; i<bps; i++) {
      if (i==bps-1 || (idp[i].first!=idp[i+1].first))
        Q[~idp[i].first] = i;
    }
    Q.back() = 0;
    this->masterQ.clear();
    this->masterQ.resize(Q.size(), 0);

    this->pCurrent = 0;
    boost::thread_group workers;
    for (int i=0; i<nthread; i++)
      workers.create_thread(boost::bind(&TopK::call64, this, bp, IdPair32(k, h), IdPair32(nthread, nPerson), boost::ref(topk_queue[i]), &Q[0]));
    workers.join_all();

    {
      for (int i=1; i<nthread; i++)
        topk_queue[0].insert(topk_queue[0].begin(), topk_queue[i].begin(), topk_queue[i].end());
      size_t qSize = std::min((size_t)k, topk_queue[0].size());
      std::partial_sort(topk_queue[0].begin(), topk_queue[0].begin()+qSize, topk_queue[0].end());
      if (qSize<topk_queue[0].size())
        topk_queue[0].resize(qSize);
    }

    std::string result_pairs = "";
    for (int i=0; i<topk_queue[0].size(); i++)
      {
        const IdPairTopK &candidate = topk_queue[0][i];
        char result_pair[64];
        sprintf(result_pair, "%llu|%llu ", candidate.second.first, candidate.second.second);
        result_pairs += result_pair;
      }
    if (!result_pairs.empty())
      result_pairs.resize(result_pairs.size()-1);
    return result_pairs;

  }

#define EXPAND(u) {                                     \
    if (!curF[u]) continue;                             \
    register BitMask bm = curF[u] & active;             \
    curF[u] = 0;                                        \
    if (!bm) continue;                                  \
    for(int k=EE[u],e=EE[u+1]; k<e; k++) {              \
      register uint32_t v = AA[k];                      \
      register BitMask nextV = bm & (~VV[v]);           \
      if (nextV) {                                      \
        nextF[v] |= nextV;                              \
        VV[v] |= nextV;                                 \
        touched |= nextV;                               \
      }                                                 \
    }                                                   \
  }

#define UPDATE() {                                      \
    active &= touched;                                  \
  }

  void call64(IdPair32 *idp, IdPair32 kh, IdPair32 stride_nPerson, TopKVector &topk_queue, const uint32_t *Q)
  {
    int maxH = kh.second;
    uint32_t personSize = this->personGraph->personSize;
    BitMask *F = (BitMask*)malloc(3*personSize*sizeof(BitMask));
    BitMask *FF[2] = {F, F + personSize};
    BitMask *VV    = F + 2*personSize;
    const uint32_t * EE = &this->personGraph->edgeIndex[0];
    const uint32_t * AA = &this->personGraph->adjacencyList[0];
    bool stop = false;
    uint32_t current_minimum = 0;
    uint32_t maxCommonTid = this->masterQ.size();
    std::vector<uint32_t> ctid(maxCommonTid);
    const uint32_t bufferSize = std::max(kh.first, 256U);
    while (true) {
      int p1Count = 0;
      uint32_t p1[NBITS];
      uint64_t person1[NBITS];
      BitMask active = 0;
      {
        this->pMutex.lock();
        {
          for (uint32_t i=current_minimum; i<maxCommonTid; i++) {
            this->masterQ[i] += ctid[i];
            if (this->masterQ[i]>=kh.first) {
              if (ctid[i]>0 && this->pCurrent>=Q[i+1])
                this->pCurrent = Q[i];
              current_minimum = i;
            }
            ctid[i] = 0;
          }
        }
        for (p1Count=0; this->pCurrent<Q[current_minimum] && p1Count<NBITS; p1Count++, this->pCurrent++) {
          p1[p1Count] = idp[this->pCurrent].second;
          active |= (1UL << p1Count);
        }
        this->pMutex.unlock();
        // fprintf(stderr, "                   \r%d/%d %u", this->pCurrent, Q[current_minimum], current_minimum);
      }
      if (active==0) break;

      memset(F, 0, 3*personSize*sizeof(BitMask));
      BitMask *curF=FF[0], *nextF=FF[0];
      for (BitMask j=0,shift=1; j<p1Count; j++, shift<<=1) {
        if (active & shift) {
          person1[j] = this->personGraph->person->denormalized(p1[j]);
          VV[p1[j]] |= shift;
          nextF[p1[j]] |= shift;
        }
      }

      uint32_t minP1 = -1;
      BitMask  touched = -1;
      // performing bfs for hop=1
      for (int hop=1; hop<2 && hop<=maxH; hop++, curF=nextF) {
        touched = 0;
        nextF = FF[hop%2];
        for (uint32_t j=0; j<p1Count; j++) {
          uint32_t u = p1[j];
          if (u<minP1) minP1 = u;
          EXPAND(u);
        }
        UPDATE();
      }
      
      // performing bfs for the rest
      for (int hop=2; hop<=maxH && touched; hop++, curF=nextF) {
        touched = 0;
        nextF = FF[hop%2];
        for (uint32_t u=0; u<personSize; u++)
          EXPAND(u);
        UPDATE();
      }

      for (uint32_t i=0, e=Q[current_minimum]; i<e; i++) {
        uint32_t p2 = idp[i].second;
        if (p2>minP1 && VV[p2]) {
          uint32_t t2End = this->personHasInterestTag->tagSize[p2];
          uint64_t person2 = this->personGraph->person->denormalized(p2);
          const Tag::IdType *tag2 = &this->personHasInterestTag->tag[this->personHasInterestTag->tagIndex[p2]];
          for (BitMask j=0,shift=1; j<p1Count; j++, shift<<=1) {
            if ((p1[j]>=p2) || (!(VV[p2]&shift))) continue;
            const Tag::IdType *tag1 = &this->personHasInterestTag->tag[this->personHasInterestTag->tagIndex[p1[j]]];
            uint32_t t1End = this->personHasInterestTag->tagSize[p1[j]];
            uint32_t common_tid = TopK::countCommon(tag1, t1End, tag2, t2End);
            if (common_tid>=current_minimum) {
              ++ctid[common_tid];
              topk_queue.push_back(IdPairTopK(~common_tid, IdPair(person1[j], person2)));
              if (topk_queue.size()>=bufferSize) {
                std::partial_sort(topk_queue.begin(), topk_queue.begin()+kh.first, topk_queue.end());
                topk_queue.resize(kh.first);
              }
            }
          }
        }
      }
      if (topk_queue.size()>=kh.first) {
        std::partial_sort(topk_queue.begin(), topk_queue.begin()+kh.first, topk_queue.end());
        topk_queue.resize(kh.first);
      }
    }
    free(F);
    // fprintf(stderr, "\n");
  }

  static inline uint32_t countCommon(const Tag::IdType *tag1, uint32_t t1End, const Tag::IdType *tag2, uint32_t t2End)
  {
    if (t1End==0 || t2End==0 || tag1[0]>tag2[t2End-1] || tag2[0]>tag1[t1End-1])
      return 0;
    register Tag::IdType t1=0, t2=0;
    register uint32_t common_tid = 0;
    while (t1<t1End && t2<t2End) {
      if (tag1[t1]<tag2[t2]) ++t1;
      else if (tag2[t2]<tag1[t1]) ++t2;
      else {
        ++common_tid;
        ++t1;
        ++t2;
      }
    }
    return common_tid;
  }
};

struct Query3Data {
  Place          * place;
  PlaceGraph     * placeGraph;
  PlaceHasPerson * placeHasPerson;
  TopK           * topK;
};

void Query3Thread::initialize()
{
  if (this->queries->empty()) return;

  Query3Data *data = new Query3Data();
  data->place = new Place(this->dataPath);
  data->placeGraph = new PlaceGraph(data->place, this->dataPath);

  const Person * person = SharedData::instance()->person;
  data->placeHasPerson = new PlaceHasPerson(person, data->place, this->dataPath);
  
  Tag *tag = SharedData::instance()->tag;
  tag->buildTagIndex();

  data->topK = new TopK(SharedData::instance()->personGraph,
                        data->placeGraph, data->placeHasPerson, tag);
  this->data = data;
}

void Query3Thread::perform()
{
  if (this->queries->empty()) return;
  this->initThread->join();
  Query3Data *data = (Query3Data*)this->data;
  for (unsigned i=0; i<this->queries->size(); i++) {
    int k, h;
    char p[128] = {};
    sscanf(this->queries->at(i)->c_str(), "query3(%d, %d, %[^)])", &k, &h, p);
    this->queries->at(i)->assign(data->topK->compute(k, h, p));
  }
  delete data->topK;
  delete data->placeHasPerson;
  delete data->placeGraph;
  delete data->place;
  delete data;
}
