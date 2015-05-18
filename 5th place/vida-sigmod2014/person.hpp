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

#ifndef PERSON_HPP
#define PERSON_HPP

#define PERSON_CSV_SORTED 1

static inline uint32_t date2uint(uint32_t year, uint32_t month, uint32_t day)
{
  return (year<<9) + (month<<5) + day;
}

struct Person
{
  typedef std::pair<uint64_t, uint32_t> IBPair;
  typedef std::vector<IBPair> IBVector;
  typedef boost::unordered_map<uint64_t, uint32_t> IdMap;
  
  IBVector info;
  IdMap    idMap;
  bool     defaultId;

  Person(const char *dataPath=0)
  {
    if (dataPath) {
      this->readPerson(dataPath);
    }
  }

  void readPerson(const char *dataPath)
  {
    this->info.clear();
    this->idMap.clear();
    std::string fn = std::string(dataPath) + "/person.csv";
    FILE *fi = fopen(fn.c_str(), "r");
    fscanf(fi, "%*[^\n]\n");
    char buffer[512];
    while (!feof(fi)) {
      IBPair ib;
      fgets(buffer, sizeof(buffer), fi);
      if (feof(fi)) break;
      sscanf(buffer, "%llu", &ib.first);
      {
        const char *p = buffer;
        int cnt = 4;
        while (cnt>0) if (*p++=='|') cnt--;
        uint32_t y, m, d;
        sscanf(p, "%u-%u-%u", &y, &m, &d);
        ib.second = date2uint(y, m, d);
      }
      this->info.push_back(ib);
    }
    fclose(fi);
    std::sort(this->info.begin(), this->info.end());
    this->defaultId = (this->info.back().first+1)==this->info.size();
    if (!defaultId)
      for (uint32_t i=0; i<this->info.size(); i++)
        this->idMap[this->info[i].first] = i;
  }
  
  inline size_t size() const
  {
    return this->info.size();
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

  inline uint64_t denormalized(uint32_t p) const
  {
    return this->info[p].first;
  }
};

#endif
