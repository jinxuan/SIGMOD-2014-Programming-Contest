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

#ifndef HSORT_HPP
#define HSORT_HPP
#include <stdint.h>

#define HSORT_SWAP(TYPE, __x, __y) { register TYPE __tmp(__x); __x = __y; __y = __tmp; }
#define HSORT_MIN_FOR_RADIX 64
namespace hsort {

  template <typename Key>
  static inline void inplaceSelectionSort(Key *A, const unsigned &n) {
    register unsigned i, j;
    for (i=0; i<n-1; i++) {
      register unsigned k(i);
      for (j=i+1; j<n; j++)
        if (A[j]<A[k])
          k = j;
      if (k!=i)
        HSORT_SWAP(Key, A[k], A[i]);
    }
  }  

  template <typename Key, typename Value>
  static inline void inplaceSelectionSortPair(void *P, const unsigned &n) {
    typedef std::pair<Key,Value> StdPair;
    StdPair *A = (StdPair*)P;
    register unsigned i, j;
    for (i=0; i<n-1; i++) {
      register unsigned k(i);
      for (j=i+1; j<n; j++)
        if (A[j]<A[k])
          k = j;
      if (k!=i)
        HSORT_SWAP(StdPair, A[k], A[i]);
    }
  }

  template <typename Key, typename Value, typename FieldType, FieldType std::pair<Key,Value>::*Field>
  static inline void inplaceSelectionSortField(void *P, const unsigned &n) {
    typedef std::pair<Key,Value> StdPair;
    StdPair *A = (StdPair*)P;
    register unsigned i, j;
    for (i=0; i<n-1; i++) {
      register unsigned k(i);
      for (j=i+1; j<n; j++)
        if (A[j].*Field<A[k].*Field)
          k = j;
      if (k!=i)
        HSORT_SWAP(StdPair, A[k], A[i]);
    }
  }

  template <typename Key, typename Value>
  static inline void inplaceSelectionSortValue(void *P, const unsigned &n) {
    typedef std::pair<Key,Value> StdPair;
    StdPair *A = (StdPair*)P;
    register unsigned i, j;
    for (i=0; i<n-1; i++) {
      register unsigned k(i);
      for (j=i+1; j<n; j++)
        if (A[j].second<A[k].second)
          k = j;
      if (k!=i)
        HSORT_SWAP(Value, A[k].second, A[i].second);
    }
  }

  template <typename Key>
  inline void inplaceRadixSort(const int &curByte, Key *A, const unsigned &n) {
    uint8_t *B = (uint8_t*)A + curByte;
    unsigned i, j, end;
    unsigned count[256] = {};

    for (i=0; i<n; i++)
      count[B[i*sizeof(Key)]]++;
    
    unsigned bucket[256];
    bucket[0] = 0;
    for (i=1; i<256; i++)
      bucket[i] = count[i-1] + bucket[i-1];

    for (i=end=0; i<256; i++) {
      if (count[i]>0) {
        end += count[i];
        for (j=bucket[i]; j<end; j++) {
          register uint8_t byte = B[j*sizeof(Key)];
          if (byte!=i) {
            register Key value(A[j]);
            do {
              register unsigned xx(bucket[byte]++);
              byte = B[xx*sizeof(Key)];
              HSORT_SWAP(Key, A[xx], value);
            } while (byte!=i);
            A[j] = value;
          }
        }
        if (curByte>0) {
          if (count[i]<HSORT_MIN_FOR_RADIX)
            hsort::inplaceSelectionSort<Key>(A+end-count[i], count[i]);
          else 
            hsort::inplaceRadixSort<Key>(curByte-1, A+end-count[i], count[i]);
        }
      }
    }
  }    

  template <typename Key, typename Value>
  inline void inplaceRadixSortKey(const int &curByte, std::pair<Key,Value> *A, const unsigned &n) {
    typedef std::pair<Key,Value> StdPair;
    uint8_t *B = (uint8_t*)A + curByte;
    unsigned i, j, end;
    unsigned count[256] = {};

    for (i=0; i<n; i++)
      count[B[i*sizeof(StdPair)]]++;
    
    unsigned bucket[256];
    bucket[0] = 0;
    for (i=1; i<256; i++)
      bucket[i] = count[i-1] + bucket[i-1];

    for (i=end=0; i<256; i++) {
      if (count[i]>0) {
        end += count[i];
        for (j=bucket[i]; j<end; j++) {
          register uint8_t byte = B[j*sizeof(StdPair)];
          if (byte!=i) {
            register StdPair value(A[j]);
            do {
              register unsigned xx(bucket[byte]++);
              byte = B[xx*sizeof(StdPair)];
              HSORT_SWAP(StdPair, A[xx], value);
            } while (byte!=i);
            A[j] = value;
          }
        }
        if (curByte>0) {
          if (count[i]<HSORT_MIN_FOR_RADIX)
            hsort::inplaceSelectionSortField<Key,Value,Key,&StdPair::first>(A+end-count[i], count[i]);
          else 
            hsort::inplaceRadixSortKey<Key,Value>(curByte-1, A+end-count[i], count[i]);
        }
      }
    }
  }  

  template <typename Key, typename Value>
  inline void inplaceRadixSortValue(const int &curByte, std::pair<Key,Value> *A, const unsigned &n) {
    typedef std::pair<Key,Value> StdPair;
    uint8_t *B = (uint8_t*)A + curByte + sizeof(Key);
    unsigned i, j, end;
    unsigned count[256] = {};

    for (i=0; i<n; i++)
      count[B[i*sizeof(StdPair)]]++;
    
    unsigned bucket[256];
    bucket[0] = 0;
    for (i=1; i<256; i++)
      bucket[i] = count[i-1] + bucket[i-1];

    for (i=end=0; i<256; i++) {
      if (count[i]>0) {
        end += count[i];
        for (j=bucket[i]; j<end; j++) {
          register uint8_t byte = B[j*sizeof(StdPair)];
          if (byte!=i) {
            register Value value(A[j].second);
            do {
              register unsigned xx(bucket[byte]++);
              byte = B[xx*sizeof(StdPair)];
              HSORT_SWAP(Value, A[xx].second, value);
            } while (byte!=i);
            A[j].second = value;
          }
        }
        if (curByte>0) {
          if (count[i]<HSORT_MIN_FOR_RADIX)
            hsort::inplaceSelectionSortValue<Key,Value>(A+end-count[i], count[i]);
          else
            hsort::inplaceRadixSortValue<Key,Value>(curByte-1, A+end-count[i], count[i]);
        }
      }
    }
  }

  template <typename Key, typename Value>
  inline void inplaceRadixSortPair(const int &curByte, std::pair<Key,Value> *A, const unsigned &n) {
    typedef std::pair<Key,Value> StdPair;
    uint8_t *B = (uint8_t*)A + curByte;
    unsigned i, j, end;
    unsigned count[256] = {};

    for (i=0; i<n; i++)
      count[B[i*sizeof(StdPair)]]++;
    
    unsigned bucket[256];
    bucket[0] = 0;
    for (i=1; i<256; i++)
      bucket[i] = count[i-1] + bucket[i-1];

    for (i=end=0; i<256; i++) {
      if (count[i]>0) {
        end += count[i];
        for (j=bucket[i]; j<end; j++) {
          register uint8_t byte = B[j*sizeof(StdPair)];
          if (byte!=i) {
            register StdPair value(A[j]);
            do {
              register unsigned xx(bucket[byte]++);
              byte = B[xx*sizeof(StdPair)];
              HSORT_SWAP(StdPair, A[xx], value);
            } while (byte!=i);
            A[j] = value;
          }
        }
        if (count[i]<HSORT_MIN_FOR_RADIX)
          hsort::inplaceSelectionSortPair<Key,Value>(A+end-count[i], count[i]);
        else {
          if (curByte>0)
            hsort::inplaceRadixSortPair<Key,Value>(curByte-1, A+end-count[i], count[i]);
          else
            hsort::inplaceRadixSortValue<Key,Value>(sizeof(Value)-1, A+end-count[i], count[i]);
        }
      }
    }
  }

  template <typename Key>
  static inline void sort(Key *A, size_t n) {
    hsort::inplaceRadixSort<Key>(sizeof(Key)-1, A, n);
  }

  template <typename Key, typename Value>
  static inline void sort(std::pair<Key, Value> *A, size_t n) {
    hsort::inplaceRadixSortPair<Key, Value>(sizeof(Key)-1, A, n);
  }

  template<class RandomAccessIterator>
  static inline void sort(RandomAccessIterator A, RandomAccessIterator B) {
    hsort::sort(&(*A), &(*B)-&(*A));
  }  

  template <typename Key, typename Value>
  static inline void sortByKey(std::pair<Key, Value> *A, size_t n) {
    hsort::inplaceRadixSortKey<Key, Value>(sizeof(Key)-1, A, n);
  }

  template<class RandomAccessIterator>
  static inline void sortByKey(RandomAccessIterator A, RandomAccessIterator B) {
    hsort::sortByKey(&(*A), &(*B)-&(*A));
  }
  
}
#undef HSORT_SWAP
#undef HSORT_MIN_FOR_RADIX

#endif
