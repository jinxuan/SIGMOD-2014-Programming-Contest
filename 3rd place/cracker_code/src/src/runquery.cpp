#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <string.h>

#define THRESHOLD 600000

int main(int argc, char** argv) {

	// get arguments for each type of query
	char* dbpath=argv[1];
	char* q1_in_path = argv[2];
	char* q1_out_path = argv[3];

	char* q2_in_path = argv[4];
	char* q2_out_path = argv[5];

	char* q3_in_path = argv[6];
	char* q3_out_path = argv[7];

	char* q4_in_path = argv[8];
	char* q4_out_path = argv[9];

	char* personcnt = argv[10];

	int pcnt = atoi(personcnt);

	pid_t query1_pid, query2_pid, query3_pid, query4_pid;


	/* pattern detection */
	char dir[100] = "";
	char str[1024];
	long delta = 10;

	strcpy(dir, dbpath);
	FILE* fp = fopen(dir, "r");

	fgets(str, 1024, fp);
	long pre, a, b, fele;
	int cnt = 10;
	char buf[10] = "10";

	while(cnt--) {
		fscanf(fp, "%ld|%ld\n", &a, &b);
		if((delta % 10) != 0){
			delta = 0;
			break;
		}
	}
	fclose(fp);
	/* end */ 


	// execute query1
	if( (query1_pid=fork()) < 0) {
		printf("1st fork error");
	} else if(query1_pid == 0) {
		if(delta == 10 )
			execl("./src/query1-seq", "", dbpath, q1_in_path, q1_out_path, personcnt, buf, (char*)0);
		else
			execl("./src/query1", "", dbpath, q1_in_path, q1_out_path, personcnt, (char*)0);
	}

	// execute query3
	if( (query3_pid=fork()) < 0) {
		printf("3rd fork error");
	} else if(query3_pid == 0) {
		if(execl("./src/query3", "", dbpath, q3_in_path, q3_out_path, personcnt, (char*)0) <0)
			printf("execute query3 failed");
	}

	// execute Query Type 1 and Query Type 2 first when graph size is bigger than THRESHOLD,
	// in order to balance CUP and IO resource and avoid memory exceeding
	if(pcnt > THRESHOLD) {
		waitpid(query3_pid, NULL, 0);
		waitpid(query1_pid, NULL, 0);
	}

	if( (query2_pid=fork()) < 0) {
		printf("2nd fork error");
	} else if(query2_pid == 0) {
		//execute query2
		if(execl("./src/query2","", dbpath, q2_in_path, q2_out_path, personcnt, (char*)0) <0)
			printf("execute query2 failed");
	}

	// execute query4
	if( (query4_pid=fork()) < 0) {
		printf("4th fork error\n");
	} else if(query4_pid == 0) {
		int ans;
		if((ans=execl("./src/query4", "", dbpath, q4_in_path, q4_out_path, personcnt, (char*)0)) <0)
			printf("execute query4 failed");

	}

	waitpid(query2_pid, NULL, 0);
	waitpid(query3_pid, NULL, 0);
	waitpid(query1_pid, NULL, 0);
	waitpid(query4_pid, NULL, 0);

	return 0;
}

