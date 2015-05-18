#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
  ios::sync_with_stdio(false);
  if (argc != 6) {
    cout << "usage: ./main query out1 out2 out3 out4" << endl;
    return 0;
  }
  ifstream qifs(argv[1]);
  ifstream oifs[4];
  
  for (int i = 0; i < 4; i++) {
    oifs[i].open(argv[i + 2]);
  }

  string line;
  while (getline(qifs, line)) {
    int q = line[5] - '1';
    getline(oifs[q], line);
    cout << line << endl;
  }
  return 0;
}
