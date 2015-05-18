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

#ifndef PERSON_GRAPH_HPP
#define PERSON_GRAPH_HPP

struct PersonGraph
{
  typedef std::pair<uint32_t, uint32_t> IdPair;
  typedef std::vector<IdPair> IdPairVector;
  typedef boost::unordered_map<IdPair, uint32_t> IdPairMap;
  typedef boost::unordered_map<uint64_t, uint32_t> U64Map;
  
  const Person *person;
  size_t        personSize;
  U64Map        relationMap;

  std::vector<uint32_t> edgeIndex;
  std::vector<uint32_t> adjacencyList;
  std::vector<uint32_t> adjacencyWeight;

  PersonGraph(const Person *P, const char *dataPath)
  {
    this->person = P;
    this->personSize = this->person->size();
    if (dataPath) {
      this->readPersonKnowsPerson(dataPath);
    }
  }

  inline uint64_t hash(const IdPair &rel) const
  {
    return rel.first*this->personSize + rel.second;
  }
  
  inline uint64_t hash(uint32_t first, uint32_t second) const
  {
    return first*this->personSize + second;
  }
  
  void readPersonKnowsPerson(const char *dataPath)
  {
    this->relationMap.clear();
    this->relationMap.rehash(this->person->size()*15/this->relationMap.max_load_factor()+1);

    IdPairVector adjacencyPair;
    {
      std::string fn = std::string(dataPath) + "/person_knows_person.csv";
      boost::iostreams::mapped_file_source file(fn);
      const char *cur = file.data();
      const char *end = cur + file.size();
      while (*cur!='\n') cur++;
      IdPair ip;
      while (++cur<end) {
        uint64_t p;
        for (p=0; *cur>='0' && *cur<='9'; cur++)
          p = p*10 + (*cur-'0');
        ip.first = this->person->normalized(p);
        for (p=0,cur++; *cur>='0' && *cur<='9'; cur++)
          p = p*10 + (*cur-'0');
        ip.second = this->person->normalized(p);     
        while (*cur!='\n') cur++;        
        adjacencyPair.push_back(ip);
      }
    }
    this->adjacencyList.resize(adjacencyPair.size());
    this->adjacencyWeight.resize(adjacencyPair.size());

    {
      hsort::sort(adjacencyPair.begin(), adjacencyPair.end());

      this->edgeIndex.clear();
      this->edgeIndex.resize(this->personSize+1, 0);

      unsigned last = 0;
      for (unsigned i=0; i<adjacencyPair.size(); i++) {
        if (i==adjacencyPair.size()-1 || adjacencyPair[i].first!=adjacencyPair[i+1].first) {
          this->edgeIndex[adjacencyPair[i].first+1] = i+1-last;
          last = i+1;
        }
        this->adjacencyList[i] = adjacencyPair[i].second;
        this->adjacencyWeight[i] = 0;
      }
      for (unsigned i=0; i<this->personSize; i++)
        this->edgeIndex[i+1] += this->edgeIndex[i];
    }
  }

  inline uint32_t normalized(IdPair rel) const
  {
    const uint32_t *alist = &this->adjacencyList[0];
    const uint32_t *begin = alist + this->edgeIndex[rel.first];
    const uint32_t *end   = alist + this->edgeIndex[rel.first+1];
    const uint32_t *v = std::lower_bound(begin, end, rel.second);
    if (v<end && *v==rel.second) {
      return v-alist;
    }
    else
      return -1;
  }
};

struct Comment
{
  typedef std::pair<uint32_t, uint32_t> IdPair;
  typedef std::pair<uint64_t, uint32_t> Id64Pair;
  typedef std::vector<uint64_t> CommentVector;
  typedef std::vector<Id64Pair> CommentPairVector;
  typedef std::vector<uint32_t> PersonVector;
  typedef std::vector<IdPair> IdPairVector;

  PersonGraph * pg;
  bool defaultId;
  bool sorted;

  static bool compareLowerBound(const Id64Pair& lhs, const uint64_t& rhs)
  {
    return lhs.first < rhs;
  }

  Comment(PersonGraph *G, const char *dataPath=0)
  {
    this->pg = G;
    if (dataPath) {
      CommentPairVector cid;
      PersonVector  pid;
      this->readComment(dataPath, pid, cid);
      this->readCommentRelation(dataPath, pid, cid);
      for (unsigned i=0; i<this->pg->personSize; i++) {
        for (unsigned k=this->pg->edgeIndex[i]; k<this->pg->edgeIndex[i+1]; k++) {
          unsigned j = this->pg->adjacencyList[k];
          if (i<j) {
            unsigned rIdx = this->pg->normalized(IdPair(j, i));
            unsigned weight = std::min(this->pg->adjacencyWeight[k],
                                       this->pg->adjacencyWeight[rIdx]);
            if (this->pg->adjacencyWeight[k]<this->pg->adjacencyWeight[rIdx])
              this->pg->adjacencyWeight[rIdx] = this->pg->adjacencyWeight[k];
            else
              this->pg->adjacencyWeight[k] = this->pg->adjacencyWeight[rIdx];
          }
        }
      }
    }
  }

  void readComment(const char *dataPath, PersonVector &pid, CommentPairVector &cid)
  {
    std::string fn = std::string(dataPath) + "/comment_hasCreator_person.csv";
    boost::iostreams::mapped_file_source file(fn);
    const char *cur = file.data();
    const char *end = cur + file.size();
    uint64_t p = 0;
    uint64_t current_c = 0;
    while (*cur!='\n') cur++;
    uint32_t lineCount = this->estimateLineCount(cur, end, file.size());
    pid.reserve(lineCount);
    this->defaultId = true;
    this->sorted = true;
    while (++cur<end) {
      if (*cur>='0' && *cur<='9')
        p = p*10 + (*cur-'0');
      else {
        if (*cur=='|') {
          if (this->defaultId && p/10!=pid.size()) {
            this->defaultId = false;
            cid.reserve(lineCount);
            for (uint64_t c=0; c<pid.size(); c++)
              cid.push_back(std::make_pair(c*10, pid[c]));
            pid.resize(0);
          }

          if (!this->defaultId) {
            if (this->sorted) {
              if (p < current_c)
                this->sorted = false;
            }
            current_c = p;
          }

        }
        else if (*cur=='\n') {
          if (this->defaultId)
            pid.push_back(this->pg->person->normalized(p));
          else
            cid.push_back(std::make_pair(current_c, this->pg->person->normalized(p)));

        }
        else
          exit(1);
        p = 0;
      }
    }
    if (!this->sorted)
      hsort::sort(cid.begin(), cid.end());
  }

  inline void resolveComment(IdPair *index, uint32_t count, const PersonVector &pid, const CommentPairVector &cid)
  {
    hsort::sortByKey(index, count);
    unsigned s=0;
    for (unsigned i=0; i<count; i++) {
      if (index[i]!=index[s]) {
        uint32_t first_pid = this->defaultId?pid[index[s].first]:cid[index[s].first].second;
        uint32_t idx = this->pg->normalized(IdPair(index[s].second, first_pid));
        if (idx!=-1)  {
          this->pg->adjacencyWeight[idx] += i-s;
        }
        s = i;
      }
    }
    uint32_t first_pid = this->defaultId?pid[index[s].first]:cid[index[s].first].second;
    uint32_t idx = this->pg->normalized(IdPair(index[s].second, first_pid));
    if (idx!=-1) 
      this->pg->adjacencyWeight[idx] += count-s;
  }

  void readCommentRelation(const char *dataPath, const PersonVector &pid, const CommentPairVector &cid)
  {
    CommentPairVector::const_iterator cb = cid.begin();
    CommentPairVector::const_iterator ce = cid.end();
    std::string fn = std::string(dataPath) + "/comment_replyOf_comment.csv";
    boost::iostreams::mapped_file_source file(fn);
    const char *cur = file.data();
    const char *end = cur + file.size();
    uint64_t p = 0;
    while (*cur!='\n') cur++;
    const int maxI = 1000000;
    std::vector<IdPair> indexV(maxI*2);
    uint32_t count = 0;
    IdPair ip;
    uint32_t *pip = (uint32_t*)&ip;
    int currentIndex = 0;
    IdPair *index[2] = {&indexV[0], &indexV[maxI]};
    boost::thread tid[2];
    while (++cur<end) {
      if (*cur>='0' && *cur<='9')
        p = p*10 + (*cur-'0');
      else {
        if (count%2==0) {
          pip[1] = this->defaultId?pid[p/10]
            :cid[std::lower_bound(cb, ce, p, Comment::compareLowerBound)-cb].second;
          ++count;
        }
        else {
          pip[0] = this->defaultId?(p/10):(std::lower_bound(cb, ce, p, Comment::compareLowerBound)-cb);
          index[currentIndex][count/2] = ip;
          if (++count>=maxI*2) {
            tid[currentIndex] = boost::thread(&Comment::resolveComment, this, index[currentIndex], count/2, boost::cref(pid), boost::cref(cid));
            currentIndex = 1-currentIndex;
            tid[currentIndex].join();
            count = 0;
          }
        }
        p = 0;
      }
    }

    if (count>0)
      this->resolveComment(index[currentIndex], count/2, boost::cref(pid), boost::cref(cid));
    tid[1-currentIndex].join();
  }
  
  static uint32_t estimateLineCount(const char *cur, const char *end, uint64_t fileSize)
  {
    uint32_t minLength=1, maxLength=1, maxCount;
    const char *i;
    for (i=cur+1; *i!='\n'; ++i) ++minLength;
    i = end-2;
    while (*i!='\n') --i;
    for (i=end-2; *i!='\n'; --i) ++maxLength;
    maxCount = fileSize*2.0/(minLength+maxLength);
    // fprintf(stderr, "Estimated line counts: minLength=%u maxLength=%u count=%u\n", minLength, maxLength, maxCount);
    return maxCount;
  }

  static void updatePersonGraph(PersonGraph *G, const char *dataPath)
  {
    Comment comm(G, dataPath);
  }

};

#endif
