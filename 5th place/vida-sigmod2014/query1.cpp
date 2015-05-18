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
#include "person_graph.hpp"

struct Query1 {
  std::vector<uint32_t> Q;
  BitSet                V;
  uint32_t            * E;
  uint32_t            * A;

  Query1(uint32_t *EE, uint32_t * AA, int maxN)
  {
    this->Q.resize(maxN);
    this->V.resize(maxN);
    this->E = EE;
    this->A = AA;
  }
  
  int hop(uint32_t src, uint32_t dst)
  {
    if (src==dst) return 0;
    int front=0, count=0;
    this->V.reset();
    this->V[src] = 1;
    this->Q[count++] = src;
    for (int hop=1; front<count; hop++) {
      for (int end=count; front<end; ++front) {
        uint32_t u = Q[front];
        for (int k=this->E[u]; k<this->E[u+1]; k++) {
          uint32_t v = this->A[k];
          if (!this->V[v]) {
            if (v==dst)
              return hop;
            this->V[v] = 1;
            this->Q[count++] = v;
          }
        }
      }
    }
    return -1;
  }

  static void runSolver(Query1Thread::QueryPair *queries, int start, int end, int stride, uint32_t *EE, uint32_t * AA, int maxN)
  {
    Query1 solver(EE, AA, maxN);
    for (int i=start; i<end; i+=stride) {
      char result[64] = {};
      sprintf(result, "%d", solver.hop(queries[i].first.second, queries[i].second.first));
      queries[i].second.second->assign(result);
    }    
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
    for (BitMask a=active; a; a&=(a-1)) {               \
      register int j = __builtin_ctzl(a);               \
      register BitMask shift = (1UL << j);              \
      if (VV[dst[j]] & shift) {                         \
        active &= ~shift;                               \
        res[j] = hop;                                   \
      }                                                 \
    }                                                   \
  }

  static void runSolver64(Query1Thread::QueryPair *queries, int start, int end, int stride, uint32_t *EE, uint32_t * AA, int maxN)
  {
    BitMask *F = (BitMask*)malloc(3*maxN*sizeof(BitMask));
    BitMask *FF[2] = {F, F + maxN};
    BitMask *VV    = F + 2*maxN;
    for (int i=start; i<end; i+=stride*NBITS) {
      memset(F, 0, 3*maxN*sizeof(BitMask));
      int      res[NBITS];
      uint32_t dst[NBITS];
      BitMask *curF=FF[0], *nextF=FF[0];
      BitMask active = -1;
      for (BitMask j=0,shift=1; j<NBITS; j++, shift<<=1) {
        uint32_t src = queries[i+j].first.second;
        dst[j] = queries[i+j].second.first;
        res[j] = -1;
        VV[src] |= shift;
        nextF[src] |= shift;
      }

      BitMask  touched = -1;
      // performing bfs for hop=1
      for (int hop=1; hop<2; hop++, curF=nextF) {
        touched = 0;
        nextF = FF[hop%2];
        for (uint32_t j=0; j<NBITS; j++) {
          uint32_t u = queries[i+j].first.second;
          EXPAND(u);
        }
        UPDATE();
      }
      
      // performing bfs for the rest
      for (int hop=2; touched; hop++, curF=nextF) {
        touched = 0;
        nextF = FF[hop%2];
        for (uint32_t u=0; u<maxN; u++)
          EXPAND(u);
        UPDATE();
      }
      
      // finalize
      for (uint32_t j=0; j<NBITS; j++) {
        char result[64] = {};
        sprintf(result, "%d", res[j]);
        queries[i+j].second.second->assign(result);
      }
    }
    free(F);
  }

  static inline void runSolverThreaded(Query1Thread::QueryPair *queries, int start, int end, uint32_t *EE, uint32_t * AA, int maxN)
  {
    int modEnd = start + (end-start)%NBITS;
    int nthread = std::min(8, modEnd-start);
    boost::thread_group workers;
    for (int i=0; i<nthread; i++)
      workers.create_thread(boost::bind(&Query1::runSolver, queries, start+i, modEnd, nthread, EE, AA, maxN));
    workers.join_all();

    nthread = std::min(8, end-modEnd);
    for (int i=0; i<nthread; i++)
      workers.create_thread(boost::bind(&Query1::runSolver64, queries, modEnd+i*NBITS, end, nthread, EE, AA, maxN));
    workers.join_all();
  }
};

void Query1Thread::initialize()
{
  if (this->queries->empty()) return;
  PersonGraph *pg = SharedData::instance()->personGraph;

  boost::thread *worker = new boost::thread(&Comment::updatePersonGraph, pg, dataPath);
  this->data = worker;
  // Comment(pg, dataPath);
  // if (pg->personSize>50000) return;

  this->sortedQueries.reserve(this->queries->size());
  for (unsigned i=0; i<this->queries->size(); i++) {
    uint64_t p1, p2;
    int x;
    sscanf(this->queries->at(i)->c_str(), "query1(%llu, %llu, %d)", &p1, &p2, &x);
    if (p1==p2)
      this->queries->at(i)->assign("0");
    else {
      this->sortedQueries.push_back(Query1Thread::QueryPair(Query1Thread::QueryKey(x, pg->person->normalized(p1)),
                                                            Query1Thread::QueryValue(pg->person->normalized(p2),
                                                                                     this->queries->at(i))));
    }
  }

  std::sort(this->sortedQueries.begin(), this->sortedQueries.end());
}

void Query1Thread::perform()
{
  if (this->queries->empty()) return;
  this->initThread->join();
  boost::thread *worker = (boost::thread*)this->data;

  PersonGraph *pg = SharedData::instance()->personGraph;
  uint32_t *EE = (uint32_t*)malloc((pg->personSize + pg->adjacencyWeight.size() + pg->adjacencyList.size())*sizeof(uint32_t));
  uint32_t *WW = EE + pg->personSize;
  uint32_t *AA = WW + pg->adjacencyWeight.size();

  uint32_t *E = &pg->edgeIndex[0];
  uint32_t *W = &pg->adjacencyWeight[0];
  uint32_t *A = &pg->adjacencyList[0];
  int prev = 0;
  int prevX = -1;
  
  int n = this->sortedQueries.size();
  for (int i=0; i<n; i++) {
    if (this->sortedQueries[i].first.first!=prevX) {
      if (prev!=i) {
        Query1::runSolverThreaded(&this->sortedQueries[0], prev, i, E, A, pg->personSize);
        prev = i;
      }
      {
        if (prevX==-1)
          worker->join();
        prevX = this->sortedQueries[i].first.first;
        int count = 0;
        uint32_t x = (prevX<0)?0:(prevX+1);
        uint32_t *prevE = E; E = EE;
        uint32_t *prevW = W; W = WW;
        uint32_t *prevA = A; A = AA;
        for (uint32_t u=0; u<pg->personSize; u++) {
          uint32_t lastEu = prevE[u];
          E[u] = count;
          for (uint32_t k=lastEu; k<prevE[u+1]; k++)
            if (prevW[k]>=x) {
              W[count] = prevW[k];
              A[count] = prevA[k];
              count++;
            }
        }
        E[pg->personSize] = count;
      }
    }
  }
  Query1::runSolverThreaded(&this->sortedQueries[0], prev, n, E, A, pg->personSize);
  free(EE);
  delete worker;
}
