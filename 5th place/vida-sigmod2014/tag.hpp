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

#ifndef TAG_HPP
#define TAG_HPP

struct Tag
{
  typedef uint32_t IdType;
  typedef std::pair<std::string, uint32_t> TagName;
  typedef std::vector<TagName> TagNameVector;
  typedef boost::unordered_map<uint64_t, uint32_t> TagIdMap;
  typedef std::pair<uint32_t, uint32_t> IdPair; 
  typedef std::vector<IdPair> IdPairVector;
  typedef std::vector<IdType> IdVector;
  typedef std::vector<uint32_t> IndexVector;
  typedef boost::unordered_map<std::string, uint32_t> TagNameMap;
 
  const Person *person;
  TagNameVector tagName;
  TagIdMap      tagMap;
  TagNameMap    tagNameMap;

  IdPairVector  byTag;
  IdPairVector  byPerson;

  IdVector      tag;
  IndexVector   tagIndex;
  IndexVector   tagSize;

  Tag(const Person *P, const char *dataPath=0)
  {
    this->person = P;
    if (dataPath) {
      this->readTag(dataPath);
      this->readTagRelation(dataPath);
    }
  }

  void readTag(const char *dataPath)
  {
    std::string fn = std::string(dataPath) + "/tag.csv";
    boost::iostreams::mapped_file_source file(fn);
    const char *cur = file.data();
    const char *end = cur + file.size();
    while (*cur++!='\n') {}
    while (cur<end) {
      {
        uint64_t p=0;
        while (*cur!='|') p = p*10 + ((*cur++)-'0');
        this->tagMap[p] = this->tagName.size();
      }
      {
        const char *start = ++cur;
        while (*cur++!='|');
        this->tagName.push_back(TagName(std::string(start, cur-start-1), tagName.size()));
      }
      while (*cur++!='\n');
    }
    this->tagNameMap.insert(this->tagName.begin(), this->tagName.end());
  }  

  void readTagRelation(const char *dataPath)
  {
    std::string fn = std::string(dataPath) + "/person_hasInterest_tag.csv";
    boost::iostreams::mapped_file_source file(fn);
    const char *cur = file.data();
    const char *end = cur + file.size();
    while (*cur!='\n') cur++;
    uint64_t p = 0;
    IdPair tp;
    this->byTag.clear();
    this->byPerson.clear();
    while (++cur<end) {
      if (*cur>='0' && *cur<='9')
        p = p*10 + (*cur-'0');
      else {
        if (*cur=='|') {
          tp.second = this->person->normalized(p);
        } else if (*cur=='\n') {
          tp.first = this->tagMap[p];
          this->byTag.push_back(tp);
          this->byPerson.push_back(IdPair(tp.second, tp.first));
        } else
          exit(1);
        p = 0;
      }
    }
  }

  void buildTagIndex()
  {
    hsort::sort(this->byPerson.begin(), this->byPerson.end());

    this->tag.resize(this->byPerson.size());
    this->tagIndex.resize(this->person->size());
    this->tagSize.resize(this->person->size());

    uint32_t current_pid = this->byPerson.front().first;
    uint32_t last_start = 0;
    for(unsigned i=0; i<this->byPerson.size(); i++) {
      if (this->byPerson[i].first != current_pid) {
        this->tagIndex[current_pid] = last_start;
        this->tagSize[current_pid] = i-last_start;
        current_pid = this->byPerson[i].first;
        last_start = i;
      }
      this->tag[i] = this->byPerson[i].second;
    }

    this->byPerson.clear();
    for (unsigned i=0; i<this->person->size(); i++) {
      if (this->tagSize[i]>0)
        this->byPerson.push_back(IdPair(~(this->tagSize[i]), i));
    }
    hsort::sort(this->byPerson.begin(), this->byPerson.end());
  }

  inline size_t size() const
  {
    return this->tagName.size();
  }

  inline uint32_t normalized(uint64_t p) const
  {
    TagIdMap::const_iterator it = this->tagMap.find(p);
    if (it!=this->tagMap.end())
      return it->second;
    return -1;
  }
};

#endif
