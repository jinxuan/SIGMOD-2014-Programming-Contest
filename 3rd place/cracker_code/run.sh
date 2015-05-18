#!/bin/bash

datapath=$1
querypath=$2

# 1. split file
./src/splitfile $querypath

personcnt=$(wc -l $datapath/person.csv | awk '{print $1}');

# 2. execute different type of query separately, parallel later
   ./src/runquery $datapath query1.in tmpquery1.out query2.in tmpquery2.out query3.in tmpquery3.out query4.in tmpquery4.out $personcnt

#    sort each type's answers to their input order
sort -n -k 1 tmpquery1.out  > query1.out
sort -n -k 1 tmpquery2.out  > query2.out
sort -n -k 1 tmpquery3.out  > query3.out
sort -n -k 1 tmpquery4.out  > query4.out

# 3. finish task by merge different type of answers
./src/mergeans
