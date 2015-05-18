#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <cstring>
#include <set>
#include <vector>
#include <queue>
#include <time.h>
#include <tr1/unordered_map>

#include <google/dense_hash_map>
#include <meminfo.hpp>
#include <sys/syscall.h>
#include <unistd.h>
#include <sched.h>

#define gettid() syscall(__NR_gettid) 
using namespace std;

using google::dense_hash_map;
using std::tr1::hash;

#define MAXLINE 50000
#define MAXTERM 100
#define PQUERY_THRESHOLD 600000

/**
 * Store an answer.
 * Contains person ID (i, j) and number of common interests.
 */
class tuple {
	public:
		int similar;
		long i, j;
		tuple():i(-1),j(-1),similar(-1) {}
		tuple(int similar_, long i_, long j_):i(i_), j(j_), similar(similar_) {}

		bool operator< (const tuple &lb) {

			if (similar != lb.similar) return similar > lb.similar;

			if (i != lb.i) return i < lb.i;

			return j < lb.j;
		}
};

/**
 * Parameters of a query.
 */
class query {
	public:
		int k;
		int h;
		string t;
		query(int k_, int h_, char* t_):k(k_), h(h_) {
			t = string(t_);
		}
};

class parameters{
	public:
		vector <int>* vp; /** valid person (in place p) */
		int batch_unit;
		int unit;
		tuple** answer_in_thread;
		bool *vpl;
		int tid;
		int h,k;
};

/**
 * Remap the person id
 */
struct eqlong {
	bool operator()(const long a, const long b) const {
		return a == b;
	}
};

inline bool operator< (const tuple &la, const tuple &lb) {

	if (la.similar != lb.similar) return la.similar > lb.similar;

	if (la.i != lb.i) return la.i < lb.i;

	return la.j < lb.j;
}

pthread_mutex_t mutex;

map <string, vector <int> > place;
tr1::unordered_map <int, int> p_ipo_p;
tr1::unordered_map <int, int> p_ili_p, o_ili_p;
tr1::unordered_map <int, vector<int> > p_wa_o, p_sa_o;
vector<int> *edge0, *p_hi_t;
int **edge;
vector <query> queries;
int queryN, curN;
int n_max;

google::dense_hash_map<long, int, hash<int>, eqlong> p2i;  /** person id --> new id */
vector <long> persons; /** new id --> person id, for finding the answer. */

/*************************************************************** Below are functions for reading the data. ***********************************************************************************/
void getplace(char *path) {
	char str[MAXLINE];
	char name[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/place.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != '|') name[i++] = *(s++);
		name[i] = 0;
		place[string(name)].push_back(strtol(id, NULL, 10));
	}
}

void getplace_isPartOf_place(char *path) {
	char str[MAXLINE];
	char id2[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/place_isPartOf_place.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != 0) id2[i++] = *(s++);
		id2[i] = 0;
		p_ipo_p[strtol(id, NULL, 10)] = strtol(id2, NULL, 10);
	}
}

void getperson_isLocatedIn_place(char *path) {
	char str[MAXLINE];
	char id2[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/person_isLocatedIn_place.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != 0) id2[i++] = *(s++);
		id2[i] = 0;
		p_ili_p[p2i[strtol(id, NULL, 10)]] = strtol(id2, NULL, 10);
	}
}

void getorganisation_isLocatedIn_place(char *path) {
	char str[MAXLINE];
	char id2[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/organisation_isLocatedIn_place.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != 0) id2[i++] = *(s++);
		id2[i] = 0;
		o_ili_p[strtol(id, NULL, 10)] = strtol(id2, NULL, 10);
	}
}

void getperson_workAt_organisation(char *path) {
	char str[MAXLINE];
	char id2[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/person_workAt_organisation.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != '|') id2[i++] = *(s++);
		id2[i] = 0;
		p_wa_o[p2i[strtol(id, NULL, 10)]].push_back(strtol(id2, NULL, 10));
	}
}

void getperson_studyAt_organisation(char *path) {
	char str[MAXLINE];
	char id2[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/person_studyAt_organisation.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != '|') id2[i++] = *(s++);
		id2[i] = 0;
		p_sa_o[p2i[strtol(id, NULL, 10)]].push_back(strtol(id2, NULL, 10));
	}
}

void getperson_hasInterest_tag(char *path) {
	char str[MAXLINE];
	char id2[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/person_hasInterest_tag.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != 0) id2[i++] = *(s++);
		id2[i] = 0;
		p_hi_t[p2i[strtol(id, NULL, 10)]].push_back(strtol(id2, NULL, 10));
	}
}

void getperson_knows_person(char *path) {
	char str[MAXLINE];
	char id2[MAXTERM];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	sprintf(file, "%s/person_knows_person.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		i = 0;
		while (*s != 0) id2[i++] = *(s++);
		id2[i] = 0;
		edge0[p2i[strtol(id, NULL, 10)]].push_back(p2i[strtol(id2, NULL, 10)]);
	}
}

void getperson(char *path) {
	char str[MAXLINE];
	char id[MAXTERM];
	char *s;

	char file[MAXTERM];
	p2i.set_empty_key(-1);
	p2i.resize(n_max);
	persons.reserve(n_max);
	sprintf(file, "%s/person.csv", path);
	freopen(file, "r", stdin);
	gets(str);
	int n = 0;
	long x;
	while (gets(str) != NULL) {
		int i = 0;
		s = str;
		while (*s != '|') id[i++] = *(s++);
		s++;
		id[i] = 0;
		x = strtol(id, NULL, 10);
		persons.push_back(x);
		p2i[x] = n++;
	}
}

/*************************************************************** End of input functions *************************************************************/

/******************************************************************** Main API **********************************************************************/

/**
 * check that whether "p" is in place "place"
 **/
inline bool checkPlace(int p, int place) {
	if (p == place) return 1;
	if (p_ipo_p.find(p) == p_ipo_p.end())
		return 0;
	else p = p_ipo_p[p];


	if (p == place) return 1;
	if (p_ipo_p.find(p) == p_ipo_p.end())
		return 0;
	else p = p_ipo_p[p];

	if (p == place) return 1;
	else return 0;
}

/**
 * Calculate the common interest tag number between two persons.
 *
 * @return Number of common interest tag.
 */
inline int calcSimilar(int i, int j) {
	int result = 0;
	int x = 0, n = p_hi_t[i].size(), m = p_hi_t[j].size();
	for (int t = 0; t < n; ++t) {
		while (x < m && p_hi_t[j][x] < p_hi_t[i][t]) x++;
		if (x < m && p_hi_t[i][t] == p_hi_t[j][x]) result++;
		if (x == m || p_hi_t[j][x] < p_hi_t[i][t]) break;
	}
	return result;
}

/**
 * Update the top-k answers.
 */
inline void calcTop(long i, long j, int similar, int k, tuple **answer) {
	tuple *tmp = new tuple(similar, i, j);
	if ((*tmp) < (*(answer[k - 1]))) {
		delete answer[k - 1];
		answer[k - 1] = tmp;
		for (int t = k - 1; t > 0; --t)
			if ((*answer[t]) < (*answer[t - 1])) {
				tmp = answer[t - 1];
				answer[t - 1] = answer[t];
				answer[t] = tmp;
			} else break;
	}
	else{
		delete tmp;
	}
}


void *batch_bfs_find(void* pnt){
	parameters* para = (parameters*)pnt;

	if(persons.size() < PQUERY_THRESHOLD){
		cpu_set_t *cpusetp;
		cpusetp = CPU_ALLOC(1);
		CPU_SET(((para->tid %4) << 1) + 1, cpusetp);
		if(sched_setaffinity(gettid(), sizeof(cpusetp), cpusetp) == -1){
		}
	}

	int start =(para-> tid) * (para->unit);
	int h = para->h;
	int k = para->k;
	int batch_unit = para->batch_unit;
	bool * vpl = para->vpl;
	vector<int>* vp = para->vp;
	tuple** answer = para->answer_in_thread;
	int up = (*vp).size();
	if(h < 1) return NULL;	

	queue <int> list;
	int n = persons.size();
	int *visitedl = new int[n + 1];
	int *dist = new int[n + 1];
	memset(visitedl, 0, sizeof(int) * (n + 1));
	int vs = 0;

	for (int idx = 0; idx < batch_unit; ++idx) {
		int i = (*vp)[start+idx];
		++vs;		
		dist[i] = 0;
		visitedl[i] = vs;
		list.push(i);
		while (!list.empty()) {
			int cur = list.front();
			list.pop();
			int size = edge[cur][0];
			for (int t = 1; t <= size; ++t) {
				int j = edge[cur][t];
				if (visitedl[j] == vs)
					continue;

				visitedl[j] = vs;
				if (dist[cur] < h - 1) {
					dist[j] = dist[cur] + 1;
					list.push(j);
				}

				if (vpl[j] && persons[i] < persons[j]) {
					if (((int)p_hi_t[i].size() >= ((*answer[k - 1])).similar) && ((int)p_hi_t[j].size() >= ((*answer[k - 1])).similar)) {
						int similar = calcSimilar(i, j);
						/*NOTE: reverse the person id here, because of the lexicographical order*/
						calcTop(persons[i], persons[j], similar, k, answer);
					}
				}
			}
		}
	}
	delete [] visitedl;
	delete [] dist;
}


/**
 * Execute a query in a single thread.
 */
void *doQuery3(void* pnt) {
	int k, h, p;
	vector <int> ps;
	int n = persons.size();
	bool* vpl = new bool[n + 1];

	while (1) {
		int current;
		pthread_mutex_lock(&mutex);
		if (curN == queryN) {
			pthread_mutex_unlock(&mutex);	
			break;
		}
		k = queries[curN].k;
		h = queries[curN].h;
		ps = place[queries[curN].t];
		current = curN++;
		pthread_mutex_unlock(&mutex);

		vector<int> vp;
		memset(vpl, 0, sizeof(bool) * (n + 1));

		for (int person = 0; person < n; ++person) {
			bool key = 0;
			for (int j = 0; j < ps.size(); ++j) {
				p = ps[j];
				if (checkPlace(p_ili_p[person], p)) {
					key = 1;
				}
				if (!key) {
					tr1::unordered_map <int, vector<int> >::iterator iter = p_wa_o.find(person); 
					if(iter != p_wa_o.end()){
						int up = (iter->second).size();
						for (int i = 0; i < up; ++i){
							if (checkPlace(o_ili_p[(iter->second)[i]], p)) {key = 1;}
						}
					}
				}
				if (!key) {
					tr1::unordered_map <int, vector<int> >::iterator iter = p_sa_o.find(person); 
					if(iter != p_sa_o.end()){
						int up = (iter->second).size();
						for(int i = 0; i < up; ++i){
							if (checkPlace(o_ili_p[(iter->second)[i]], p)) {key = 1;}
						}
					}
				}
			}
			if (!key) continue;
			vp.push_back(person);
			vpl[person] = 1;
		}

		if(n < PQUERY_THRESHOLD) {
			/* Use bfs to find every person's h-neighbor. */
			tuple *answer[k + 1];
			for (int i = 0; i < k; ++i)
				answer[i] = new tuple(-1, -1, -1);

			if (h > 0){ 
				queue <int> list;
				int *visitedl = new int[n + 1];
				int *dist = new int[n + 1];
				memset(visitedl, 0, sizeof(int) * (n + 1));
				int vs = 0;

				for (int idx = 0; idx < vp.size(); ++idx) {
					int i = vp[idx];
					++vs;		
					dist[i] = 0;
					visitedl[i] = vs;
					list.push(i);
					while (!list.empty()) {
						int cur = list.front();
						list.pop();
						int size = edge[cur][0];
						for (int t = 1; t <= size; ++t) {
							int j = edge[cur][t];
							if (visitedl[j] == vs)
								continue;

							visitedl[j] = vs;
							if (dist[cur] < h - 1) {
								dist[j] = dist[cur] + 1;
								list.push(j);
							}

							/* Get the number of common interest tag and update the top-k answer. */
							if (vpl[j] && persons[i] < persons[j]) {
								if (((int)p_hi_t[i].size() >= ((*answer[k - 1])).similar) && ((int)p_hi_t[j].size() >= ((*answer[k - 1])).similar)) {
									int similar = calcSimilar(i, j);
									calcTop(persons[i], persons[j], similar, k, answer);
								}
							}
						}
					}
				}
				delete[] visitedl;
				delete[] dist;
			}	

			char result[MAXLINE];
			char tmp[MAXTERM];
			sprintf(result, "%d", current);
			tuple *iter;
			for (int i = 0; i < k; ++i) {
				iter = answer[i];
				if ((*iter).similar < 0) break;
				sprintf(tmp, " %ld|%ld", (*iter).i, (*iter).j);
				strcat(result, tmp);
				delete answer[i];
			}
			printf("%s\n", result);
		}
		else {
			/* Use thread parallelism when graph size is too big. */
			int thread_num = 8;
			int unit = (vp.size() + thread_num - 1) / thread_num;
			pthread_t pt[thread_num];
			parameters* th_para = new parameters[thread_num];

			for (int i = 0; i < thread_num; ++i) {
				th_para[i].vp = &vp;
				th_para[i].unit = unit;
				th_para[i].batch_unit = (i+1)*unit > vp.size() ? (vp.size() - i*unit) : unit;
				th_para[i].tid = i;
				th_para[i].answer_in_thread = new tuple*[k+1];
				for (int j = 0; j < k; ++j)
					th_para[i].answer_in_thread[j] = new tuple();
				th_para[i].vpl = vpl;
				th_para[i].h = h;
				th_para[i].k = k;

				pthread_create(&pt[i], NULL, batch_bfs_find, (void*)(th_para+i));
			}

			for (int i = 0; i < thread_num; ++i) {
				pthread_join(pt[i], NULL);
			}

			/* aggregate answer */
			for(int i = 1; i < thread_num; ++i){
				tuple** iter = th_para[i].answer_in_thread;
				for(int j = 0; j < k; ++j){
					calcTop(iter[j]->i, iter[j]->j, iter[j]->similar, k, th_para[0].answer_in_thread);
					delete iter[j];
				}
				delete [] th_para[i].answer_in_thread;
			}

			char result[MAXLINE];
			char tmp[MAXTERM];
			sprintf(result, "%d", current);
			tuple *iter;
			for (int i = 0; i < k; ++i) {
				iter = *(th_para[0].answer_in_thread + i);
				if ((*iter).similar < 0) break;
				sprintf(tmp, " %ld|%ld", (*iter).i, (*iter).j);
				strcat(result, tmp);
				delete iter;
			}
			printf("%s\n", result);
			delete[] th_para;   
		}
		vp.clear();
	}

	delete[] vpl;
}

int main(int argc, char** argv) {
	char *path = argv[1];
	int n_max = atoi(argv[4]);
	getperson(path);

	int n = persons.size();
	edge0 = new vector <int>[n + 1]; 
	p_hi_t = new vector <int>[n + 1]; 

	/* Initialization */
	getplace(path);
	getplace_isPartOf_place(path);
	getperson_isLocatedIn_place(path);
	getorganisation_isLocatedIn_place(path);
	getperson_workAt_organisation(path);
	getperson_studyAt_organisation(path);
	getperson_hasInterest_tag(path);
	getperson_knows_person(path);

	edge = new int*[n + 1];
	for (int i = 0; i < n; ++i) {
		int size = edge0[i].size();
		edge[i] = new int[size + 2];
		edge[i][0] = size;
		memmove(&edge[i][1], edge0[i].data(), sizeof(int) * (size));
	}
	delete[] edge0;

	for (int i = 0; i < n; ++i) {
		sort(p_hi_t[i].begin(), p_hi_t[i].end());
	}

	freopen(argv[2], "r", stdin);
	freopen(argv[3], "w", stdout);
	char str[MAXLINE];
	void *status;


	while (gets(str) != NULL) {
		char *s = strchr(str, '(');
		s++;
		int n = strlen(s);
		for (int i = 0; i < n; ++i)
			if (s[i] == ',') s[i] = ' ';

		s[n - 1] = 0;
		int k, h;
		char *t = new char[MAXLINE];
		sscanf(s, "%d %d %s", &k, &h, t);
		queryN++;
		queries.push_back(query(k, h, t));
	}
	int THREADN = 10;

	if(persons.size() > PQUERY_THRESHOLD)
		THREADN  = 8;

	pthread_t pt[THREADN];
	pthread_mutex_init(&mutex, NULL);	

	for (int i = 0; i < THREADN; ++i) {
		/* execute query here */
		pthread_create(&pt[i], NULL, doQuery3, new int(i));
	}

	for (int i = 0; i < THREADN; ++i) {
		pthread_join(pt[i], NULL);
	}

	n = persons.size();
	for (int i = 0; i < n; ++i)
		delete[] edge[i];
	delete[] edge;
}
