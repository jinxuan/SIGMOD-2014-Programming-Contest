#include <cstdio>
#include <cstdlib>

using namespace std;

#define MAX_LEN 50000


char metafile[32] = "query.split.meta";
char q1outfile[15] = "query1.out";
char q2outfile[15] = "query2.out";
char q3outfile[15] = "query3.out";
char q4outfile[15] = "query4.out";

int main(int argc, char** args) {
	FILE* inmeta = fopen(metafile, "rb");

	FILE* ansq[5];
	ansq[1] = fopen(q1outfile,"r");
	ansq[2] = fopen(q2outfile,"r");
	ansq[3] = fopen(q3outfile,"r");
	ansq[4] = fopen(q4outfile,"r");

	int qt, cnt;
	char ansline[MAX_LEN];
	char*pt;

	int qid;

	while(fscanf(inmeta, "%d %d\n", &qt, &cnt) != EOF) {
		for(int i = 0; i < cnt; ++i) {
			fgets(ansline, MAX_LEN, ansq[qt]);
			pt=ansline;
			while((*pt) != ' ' && (*pt) != '\n') pt++;
			if((*pt) == ' ')
				printf("%s", pt+1);
			else
				printf("\n");
		}
	}
	return 0;
}
