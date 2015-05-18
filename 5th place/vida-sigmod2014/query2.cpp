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

struct Query2 {
  typedef std::pair<uint32_t, uint32_t> IdPair;

  Tag         * tag;

  std::vector<uint32_t> resultIndex;
  uint32_t            * resultCount;
  

  Query2(Tag *T): tag(T)
  {
    this->computeAll();
  }
  
  void computeAll()
  { 
    std::sort(this->tag->tagName.begin(), this->tag->tagName.end(), std::greater<Tag::TagName>());
    hsort::sort(this->tag->byTag.begin(), this->tag->byTag.end());
    
    PersonGraph   * pg = SharedData::instance()->personGraph;
    
    if (this->tag->byTag.size()>0) {
      this->tag->byTag.push_back(IdPair(-1,-1));
        
      this->resultIndex.clear();
      this->resultIndex.resize(this->tag->tagName.size()*2, 0);
      this->resultCount = &this->resultIndex[this->tag->tagName.size()];

      unsigned start = 0;
      uint32_t tid = this->tag->byTag[start].first;
      boost::unordered_set<int32_t> vSet;
      std::vector<uint32_t> cc(pg->personSize, -1);
      std::vector<uint32_t> next(pg->personSize, 0);
      
      for (unsigned i=0; i<this->tag->byTag.size(); i++) {
        if (this->tag->byTag[i].first!=tid) {

          for (unsigned j=start; j<i; j++)
            this->tag->byTag[j].first = pg->person->info[this->tag->byTag[j].second].second;

          std::sort(&this->tag->byTag[start], &this->tag->byTag[i], std::greater<IdPair>());

          vSet.clear();
          vSet.rehash((i-start)/vSet.max_load_factor()+1);
          this->resultIndex[tid] = start;
          this->resultCount[tid] = i-start;
          int maxCnt = 0;
          for (unsigned j=start; j<i; j++) {
            uint32_t u = this->tag->byTag[j].second;
            cc[u] = u;
            next[u] = -1;
            if (maxCnt<1) maxCnt = 1;

            for (uint32_t k=pg->edgeIndex[u]; k<pg->edgeIndex[u+1]; k++) {
              uint32_t v = pg->adjacencyList[k];
              if (vSet.find(v)!=vSet.end() && cc[u]!=cc[v]) {
                uint32_t n;
                uint32_t cnt=1;
                for (n=cc[v]; next[n]!=-1; n=next[n],cnt++);
                next[n] = cc[u];
                for (n=cc[u]; n!=-1; n=next[n],cnt++)
                  cc[n] = cc[v];
                if (maxCnt<cnt) maxCnt = cnt;
              }
            }
            vSet.insert(u);
            this->tag->byTag[j].second = maxCnt;
          }

          std::reverse(&this->tag->byTag[start], &this->tag->byTag[i]);

          start = i;
          tid = this->tag->byTag[start].first;
        }
      }

      this->tag->byTag.pop_back();
    }
  }

  void solve(std::string& query)
  {
    uint32_t k, year, month, day, bd;
    sscanf(query.c_str(), "query2(%u, %u-%u-%u)\n", &k, &year, &month, &day);
    bd = date2uint(year, month, day);
    query.clear();
    std::vector<IdPair> candidates(this->tag->tagName.size());
    for (unsigned i=0; i<this->tag->tagName.size(); i++) {
      uint32_t cnt = 0;
      int j = this->tag->tagName[i].second;
      const IdPair *p = &this->tag->byTag[this->resultIndex[j]];
      if (this->resultCount[j]>0 && p[this->resultCount[j]-1].first>=bd) {
        cnt = std::lower_bound(p, p+this->resultCount[j], IdPair(bd, 0))->second;
      }
      candidates[i] = IdPair(cnt, i);
    }

    k = std::min(k,(uint32_t)candidates.size());
    std::partial_sort(candidates.begin(), candidates.begin()+k, candidates.end(), std::greater<IdPair>());
    candidates.resize(k);

    for (unsigned i=0; i<k; i++) {
      query += this->tag->tagName[candidates[i].second].first;
      if (k>1) query += " ";
    }
  }

};

void Query2Thread::initialize()
{
  if (this->queries->empty()) return;
  Query2 solver(SharedData::instance()->tag);
  for (unsigned i=0; i<this->queries->size(); i++) {
    solver.solve(*this->queries->at(i));
  }
}

void Query2Thread::perform()
{
  if (this->queries->empty()) return;
  this->initThread->join();
}
