#include <cstdio>
#include <cstdlib>

using namespace std;

#define MAX_LEN 256

char metafile[32] = "query.split.meta";

char query1file[10] = "query1.in";
char query2file[10] = "query2.in";
char query3file[10] = "query3.in";
char query4file[10] = "query4.in";

/**
 * Split query file by different type of queries.
 * Each type of queries will be handled by a single process. 
 */
int main(int argc, char** args) {
	char* querypath = args[1];
	FILE* infp = fopen(querypath, "r");
	FILE* outmeta = fopen(metafile, "wb");

	FILE* outq[5];
	outq[1] = fopen(query1file,"w");
	outq[2] = fopen(query2file,"w");
	outq[3] = fopen(query3file,"w");
	outq[4] = fopen(query4file,"w");

	char inputline[MAX_LEN];
	int preqt = 0;
	int curqt = 0;
	int preline = -1;
	int curline = -1;

	if(fgets(inputline, MAX_LEN, infp) != NULL) {
		preqt = inputline[5] - 48; //char to int
		preline = 0;
		curline = 0;
		fputs(inputline, outq[preqt]);
	}

	while(fgets(inputline, MAX_LEN, infp) != NULL) {
		curline++;
		curqt = inputline[5] - 48; //char to int

		fputs(inputline, outq[curqt]);

		if(curqt != preqt){
			fprintf(outmeta, "%d %d\n",preqt,curline-preline); //querytype, number
			preqt = curqt;
			preline=curline;
		}
	}

	fprintf(outmeta, "%d %d\n",preqt,curline-preline+1); //querytype, number

	fclose(infp);
	fclose(outq[1]);
	fclose(outq[2]);
	fclose(outq[3]);
	fclose(outq[4]);
	fclose(outmeta);
	return 0;
}
