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

#ifndef QUERY_HPP
#define QUERY_HPP

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include <boost/timer/timer.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "hsort.hpp"
#include "person.hpp"
#include "person_graph.hpp"
#include "tag.hpp"

#define NBITS 64
typedef uint64_t BitMask;
typedef boost::dynamic_bitset<uint64_t> BitSet;
typedef std::pair<int, std::string> Query;
typedef std::vector<std::string> QueryVec;
typedef std::pair<std::string, int> QueryPair;
typedef std::vector<QueryPair> QueryPairVec;
typedef std::vector<std::string*> QueryVecPtr;
typedef std::pair<std::string*, int> QueryPairPtr;
typedef std::vector<QueryPairPtr> QueryPairVecPtr;

struct SharedData {
  SharedData();
  void initialize(const char *path);
  static SharedData * instance();

  const char  * dataPath;
  Person      * person;
  PersonGraph * personGraph;
  Tag         * tag;
};

struct QueryThread {
  QueryThread(QueryVecPtr *Q): dataPath(SharedData::instance()->dataPath), data(0), queries(Q), initThread(0) {}
  virtual ~QueryThread();
  
  void startInitialization();
  void startPerforming();

  virtual void initialize() = 0;
  virtual void perform() = 0;

  const char    * dataPath;
  void          * data;
  QueryVecPtr   * queries;
  boost::thread * initThread;
};

struct Query1Thread : public QueryThread {
  Query1Thread(QueryVecPtr *Q): QueryThread(Q) {}
  void initialize();
  void perform();
  
  typedef std::pair<int, uint32_t> QueryKey;
  typedef std::pair<uint32_t, std::string*> QueryValue;
  typedef std::pair<QueryKey, QueryValue> QueryPair;
  typedef std::vector<QueryPair> QueryPairVector;
  QueryPairVector sortedQueries;
};

struct Query2Thread : public QueryThread {
  Query2Thread(QueryVecPtr *Q): QueryThread(Q) {}
  void initialize();
  void perform();
};

struct Query3Thread : public QueryThread {
  Query3Thread(QueryVecPtr *Q): QueryThread(Q) {}
  void initialize();
  void perform();
};

struct Query4Thread : public QueryThread {
  Query4Thread(QueryVecPtr *Q): QueryThread(Q) {}
  void initialize();
  void perform();
};

#endif
