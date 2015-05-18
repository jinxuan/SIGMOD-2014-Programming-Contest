CXX = g++ -std=c++0x -Wall -Wextra -Wunused-function
CXXARCHFLAGS = -march=core2 -mcx16 -msahf -msse4.1 --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=6144 -mtune=core2
#CXXFLAGS = -O3 -fopenmp ${CXXARCHFLAGS}
CXXFLAGS = -O3  -fopenmp ${CXXARCHFLAGS}
all: query1 query2 query3 query4 main

main: main.cpp
	${CXX} ${CXXFLAGS} -o $@ $^

query1: query1.cc
	${CXX} ${CXXFLAGS} -o $@ $^

query2: query2.cc
	${CXX} ${CXXFLAGS} -o $@ $^

query3: query3.cc
	${CXX} ${CXXFLAGS} -o $@ $^

query4: query4.cc
	${CXX} ${CXXFLAGS} -o $@ $^



small: query1 query3
	echo "****************************** Query1 ******************************"
	./query1 outputDir-1k/ data/1k-queries.txt > data/1k-Q1-outputs.txt
	diff -sq --strip-trailing-cr data/1k-Q1-outputs.txt data/1k-Q1-answers.txt
	echo "****************************** Query2 ******************************"
	./query2 outputDir-1k/ data/1k-queries.txt > data/1k-Q2-outputs.txt
	diff -sq --strip-trailing-cr data/1k-Q2-outputs.txt data/1k-Q2-answers.txt
	echo "****************************** Query3 ******************************"
	./query3 outputDir-1k/ data/1k-queries.txt > data/1k-Q3-outputs.txt
	diff -sq --strip-trailing-cr data/1k-Q3-outputs.txt data/1k-Q3-answers.txt
	echo "****************************** Query4 ******************************"
	./query4 outputDir-1k/ data/1k-queries.txt > data/1k-Q4-outputs.txt
	diff -sq --strip-trailing-cr data/1k-Q4-outputs.txt data/1k-Q4-answers.txt

medium: query1 query3
	echo "****************************** Query1 ******************************"
	./query1 outputDir-10k/ data/10k-queries.txt > data/10k-Q1-outputs.txt
	diff -sq --strip-trailing-cr data/10k-Q1-outputs.txt data/10k-Q1-answers.txt
	echo "****************************** Query2 ******************************"
	./query2 outputDir-10k/ data/10k-queries.txt > data/10k-Q2-outputs.txt
	diff -sq --strip-trailing-cr data/10k-Q2-outputs.txt data/10k-Q2-answers.txt
	echo "****************************** Query3 ******************************"
	./query3 outputDir-10k/ data/10k-queries.txt > data/10k-Q3-outputs.txt
	diff -sq --strip-trailing-cr data/10k-Q3-outputs.txt data/10k-Q3-answers.txt
	echo "****************************** Query4 ******************************"
	./query4 outputDir-10k/ data/10k-queries.txt > data/10k-Q4-outputs.txt
	diff -sq --strip-trailing-cr data/10k-Q4-outputs.txt data/10k-Q4-answers.txt

.PHONY: clean

clean:
	rm -rf query1 query2 query3 query4 main
