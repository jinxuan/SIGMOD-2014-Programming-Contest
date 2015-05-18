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
#include <boost/thread.hpp>
#include "person_graph.hpp"

#define DEBUG 0

SharedData * SharedData::instance()
{
  static SharedData *global = 0;
  if (!global)
    global = new SharedData;
  return global;
}

SharedData::SharedData(): person(0), personGraph(0), dataPath(0)
{
}

void SharedData::initialize(const char *path)
{
  this->dataPath = path;
  this->person = new Person(this->dataPath);
  this->tag = new Tag(this->person, this->dataPath);
  this->personGraph = new PersonGraph(this->person, this->dataPath);
}

QueryThread::~QueryThread()
{
  delete this->initThread;
}

void QueryThread::startInitialization()
{
  this->initThread = new boost::thread(&QueryThread::initialize, this);
}

void QueryThread::startPerforming()
{
  this->initThread->join();
  this->perform();
}

bool compareQueryResults(const QueryPairPtr& lhs, const QueryPairPtr& rhs)
{
  return lhs.second < rhs.second;
}

int main(int argc, char**argv)
{
  const char *dataPath  = "../data/outputDir-1k";
  const char *queryFile = "../queries/1k-queries.txt";

  // Read general arguments
  for (int i=1; i<argc; i++) {
    if (strcmp(argv[i], "--help")==0) {
      fprintf(stderr, "Usage: %s [options]\n", argv[0]);
      fprintf(stderr, "Options:\n");
      fprintf(stderr, "   --help\n");
      fprintf(stderr, "   --data     <data path>   [default: %s]\n", dataPath);
      fprintf(stderr, "   --query    <query file>  [default: %s]\n", queryFile);
      exit(0);
    }
    else if (strcmp(argv[i], "--data")==0) {
      dataPath = argv[++i];
    }
    else if (strcmp(argv[i], "--query")==0) {
      queryFile = argv[++i];
    }
    else {
      fprintf(stderr, "Ignored unknown arguments '%s'.\n", argv[i++]);
    }
  }

  QueryPairVec queries;
  
  FILE *fi = fopen(queryFile, "r");
  int counter = 0;
  while (!feof(fi)) {
    char query[256] = {};
    fgets(query, sizeof(query), fi);
    std::string queryStr = std::string(query);
    if (queryStr.length() < 5)
      continue;

    queries.push_back(std::make_pair(queryStr, counter++));
    if (queryStr[5]<'1' || queryStr[5]>'4') {
      std::cerr << "Sorry, wrong query: " << queryStr << std::endl;
      exit(1);
    }
  }
  fclose(fi);

  std::sort(queries.begin(), queries.end());

  QueryVecPtr pQueries[4];
  QueryPairVecPtr queryResults(queries.size());
  std::string* previousQuery = &(queries[0].first);
  pQueries[queries[0].first[5]-'1'].push_back(&(queries[0].first));
  queryResults[0] = std::make_pair(&(queries[0].first), queries[0].second);
  for (unsigned i=1; i<queries.size(); i++)
    {
      if (queries[i].first.compare(*previousQuery) == 0)
        queryResults[i] = std::make_pair(previousQuery, queries[i].second);
      else
        {
          pQueries[queries[i].first[5]-'1'].push_back(&(queries[i].first));
          queryResults[i] = std::make_pair(&(queries[i].first), queries[i].second);
          previousQuery = &(queries[i].first);
        }
    }

  {
    boost::timer::cpu_timer timer;
    SharedData *data = SharedData::instance();
    data->initialize(dataPath);

#if DEBUG
    if (data->person->size()>50000)
      fprintf(stdout, "0:%d ", (int)(timer.elapsed().wall*2e-6));
#endif

    QueryThread *query[5] = { 0,
                              new Query1Thread(&pQueries[0]),
                              new Query2Thread(&pQueries[1]),
                              new Query3Thread(&pQueries[2]),
                              new Query4Thread(&pQueries[3])};

    query[3]->startInitialization();
    query[4]->startInitialization();
    query[2]->startInitialization();
    query[1]->startInitialization();
#if DEBUG
    query[3]->initThread->join();
    if (data->person->size()>50000)
      fprintf(stdout, "3*:%d ", (int)(timer.elapsed().wall*2e-6));
#endif

    query[3]->perform();
    
#if DEBUG    
    if (data->person->size()>50000)
      fprintf(stdout, "3:%d ", (int)(timer.elapsed().wall*2e-6));    
    query[4]->initThread->join();
    if (data->person->size()>50000)
      fprintf(stdout, "4*:%d ", (int)(timer.elapsed().wall*2e-6));
#endif

    query[4]->perform();

#if DEBUG    
    if (data->person->size()>50000)
      fprintf(stdout, "4:%d ", (int)(timer.elapsed().wall*2e-6));
#endif

    query[2]->perform();

#if DEBUG
    if (data->person->size()>50000)
      fprintf(stdout, "2:%d ", (int)(timer.elapsed().wall*2e-6));

    query[1]->initThread->join();

    if (data->person->size()>50000) {
      fprintf(stdout, "1*:%d\n", (int)(timer.elapsed().wall*2e-6));
    }
#endif
    
    query[1]->perform();

#if DEBUG
    if (data->person->size()>50000) {
      fprintf(stdout, "1:%d\n", (int)(timer.elapsed().wall*2e-6));
      return 0;
    }
#endif
  }

  std::sort(queryResults.begin(), queryResults.end(), compareQueryResults);

  for (unsigned i=0; i<queryResults.size(); i++)
    {
      boost::trim_right(*(queryResults[i].first));
      std::cout << *(queryResults[i].first) << std::endl;
    }

  return 0; 
}
