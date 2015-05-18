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
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/wait.h>

#include <google/dense_hash_map>

using namespace std;

using google::dense_hash_map;
using std::tr1::hash;

#include <io_mmap.hpp>

#define MAXLINE 50000
#define MAXTERM 100
#define THREADN 10
#define MAXDISTANCE (5)
#define phi (0.77351)
#define NMAP 128

#define EPSILON 1e-6

/******************************************************************************* Class and struct deffinition **********************************************************************/

/**
 * record the closeness cetrality of vertex id
 */
class vertex_centrality {
	public:
		long id;
		double closeness;
		vertex_centrality() { }
		vertex_centrality(double closeness_, int id_):id(id_), closeness(closeness_) { }

};

/**
 * record the average distance of vertex id
 */
class vertex_avgdist {
	public:
		int id;
		double avgdist;
		vertex_avgdist() { }
		vertex_avgdist(double avgdist_, int id_):id(id_), avgdist(avgdist_) { }

};

/**
 * record the vertex id and its degree
 */
class vertex_degree {
	public:
		int id;
		int degree;
		vertex_degree() { }
		vertex_degree(int degree_, int id_):id(id_), degree(degree_) { }

};

/**
 * query and database related structure
 */
class query {
	public:
		int k;
		int t;
		query(int k_, int s_):k(k_), t(s_) {
		}
};

/**
 * remap the person id
 */
struct eqlong {
	bool operator()(const long a, const long b) const {
		return a == b;
	}
};

/************************************************************************* Override operator ***************************************************************************/

double abs(double t) {
	if (t < 0) return -t;
	return t;
}

inline bool operator< (const vertex_centrality &la, const vertex_centrality &lb) {
	if (abs(la.closeness - lb.closeness) > EPSILON) 
		return la.closeness > lb.closeness;
	return la.id < lb.id;

}

inline bool operator< (const vertex_avgdist &la, const vertex_avgdist &lb) {
	if (abs(la.avgdist - lb.avgdist) > EPSILON) 
		return la.avgdist < lb.avgdist;
	return la.id < lb.id;

}

inline bool operator< (const vertex_degree &la, const vertex_degree &lb) {
	if (la.degree != lb.degree) 
		return la.degree < lb.degree;
	return la.id < lb.id;

}

/************************************************************************* Variables ***************************************************************************/

char dir[MAXTERM]; /** Initial IO related */
google::dense_hash_map<long, int, hash<int>, eqlong> p2i;   //person id --> new id
vector <long> persons; // new id --> person id, for finding the answer.
int n_max;

vector<int>* edge;
map <int, set<int> > f_ht_t;
vector<int>* p_hf_f;
map <string, int> tag;
vector <query> queries;
int queryN, curN;
time_t start_time;

/* thread mutex */
pthread_mutex_t mutex;

/* declaration */
void *doQuery4(void *pnt);
void runQuery(int qid, int k, int t);
void topkrank(int qid, int k, vector<int> valid_person_list, vector<int> *person_graph, int num_person, int total_num_person);
void cyclic_calculation(int qid, set<vertex_centrality>& answer, int k, vector<vertex_degree>& cc_vertices, int cc_size, 
		vector <int> * person_graph, int num_person, int total_num_person);
void update_answer(set<vertex_centrality>& answer, int k, double centality, long vid);
void bfs(int start, int* estimated_distances, vector <int> * person_graph, int total_num_person, int& estimated_radius);


/************************************************************************* Several hash functions ***************************************************************************/
inline unsigned int hash1(unsigned int a) {
	a = (a + 0x7ed55d16) + (a << 12);
	a = (a ^ 0xc761c23c) ^ (a >> 19);
	a = (a + 0x165667b1) + (a << 5);
	a = (a + 0xd3a2646c) ^ (a << 9);
	a = (a + 0xfd7046c5) + (a << 3);
	a = (a ^ 0xb55a4f09) ^ (a >> 16);
	return a;
}

inline uint32_t hash2( uint32_t a) {
	a += ~(a<<15);
	a ^=  (a>>10);
	a +=  (a<<3);
	a ^=  (a>>6);
	a += ~(a<<11);
	a ^=  (a>>16);
}


inline uint32_t popcnt( uint32_t x ) {
	x -= ((x >> 1) & 0x55555555);
	x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
	x = (((x >> 4) + x) & 0x0f0f0f0f);
	x += (x >> 8);
	x += (x >> 16);
	return x & 0x0000003f;
}

uint32_t hash3( uint32_t a) {
	a = (a^0xdeadbeef) + (a<<4);
	a = a ^ (a>>10);
	a = a + (a<<7);
	a = a ^ (a>>13);
	return a;
}

/************************************************************************* Functions for reading the data ***************************************************************************/
void getperson(char *dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;
	p2i.set_empty_key(-1);
	p2i.resize(n_max);
	persons.reserve(n_max);

	if(mmap_file_reader.open_mmap(strcat(dir,"/person.csv")) == false)
		return;
	mmap_file_reader.skipfirstline();

	long id;
	int n = 0;
	while (mmap_file_reader.hasNext()) {
		id = mmap_file_reader.getLong();
		persons.push_back(id);
		p2i[id] = n++;
		mmap_file_reader.skipline();
	}
	mmap_file_reader.close_mmap();
}

void getperson_knows_person(char *dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;
	if(mmap_file_reader.open_mmap(strcat(dir,"/person_knows_person.csv")) == false)
		return;
	mmap_file_reader.skipfirstline();

	long id, id2;
	edge = new vector<int>[persons.size()+1];
	int n = persons.size();
	for (int i = 0; i < n; ++i)
		edge[i].reserve(100);
	while(mmap_file_reader.hasNext()) {
		id = mmap_file_reader.getLong();
		id2 = mmap_file_reader.getLong();
		edge[p2i[id]].push_back(p2i[id2]);
	}
	mmap_file_reader.close_mmap();
}

void getforum_hasTag_tag(char *dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;
	if(mmap_file_reader.open_mmap(strcat(dir,"/forum_hasTag_tag.csv")) == false)
		return;
	mmap_file_reader.skipfirstline();

	int id, id2;
	while (mmap_file_reader.hasNext()) {
		id = mmap_file_reader.getInt();
		id2 = mmap_file_reader.getInt();
		f_ht_t[id].insert(id2);
	}
	mmap_file_reader.close_mmap();
}

void getforum_hasMember_person(char *dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;

	if(mmap_file_reader.open_mmap(strcat(dir,"/forum_hasMember_person.csv")) == false)
		return;

	long limit = 5L * 1024 * 1024 * 1024; //5G
	long offset = 0, readbytes = 0;

	readbytes = mmap_file_reader.skipfirstline();

	int id;
	long id2;

	p_hf_f = new vector<int>[persons.size()+1];
	int n = persons.size();
	for (int i = 0; i < n; ++i)
		p_hf_f[i].reserve(100);

	while (mmap_file_reader.hasNext()) {
		if(readbytes > limit) {
			offset += readbytes;
			readbytes = 0;
			bool flag = mmap_file_reader.remap(dir, offset);
		}

		readbytes += mmap_file_reader.getInt(id);
		readbytes += mmap_file_reader.getLong(id2);

		p_hf_f[p2i[id2]].push_back(id);
		readbytes += mmap_file_reader.skipline();
	}
	mmap_file_reader.close_mmap();
}

void gettag(char* dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;
	if(mmap_file_reader.open_mmap(strcat(dir,"/tag.csv")) == false)
		return;
	mmap_file_reader.skipfirstline();

	int id;
	char name[MAXTERM];

	while (mmap_file_reader.hasNext()) {
		id = mmap_file_reader.getInt();
		mmap_file_reader.getString(name);
		tag[string(name)] = id;
		mmap_file_reader.skipline();
	}
	mmap_file_reader.close_mmap();
}


/************************************************************************* Main API ***************************************************************************/
/**
 * Main body
 */
int main(int argc, char** argv) {
	/* initialize database */
	start_time = time(NULL);
	char *path = argv[1];
	n_max = atoi(argv[4]);

	getperson(path);
	getperson_knows_person(path);
	getforum_hasTag_tag(path);
	getforum_hasMember_person(path);
	gettag(path);

	pid_t query1_pid;
	if(persons.size() > 6000000) {
		int delta = atoi(argv[5]);
		if( (query1_pid=fork()) < 0) {
			printf("1st fork error");
		} else if(query1_pid == 0) {
			/* execute query1 */
			if(delta == 10) execl("./src/query1-seq", "", path, argv[6], argv[7], argv[4], argv[5], (char*)0);
				else execl("./src/query1", "", path, argv[6], argv[7], argv[4], (char*)0);
		}
	}


	/* read the queries */
	freopen(argv[2], "r", stdin);
	freopen(argv[3], "w", stdout);
	FILE* fp = fopen(argv[3], "w");
	char str[MAXLINE];
	void *status;
	char t[MAXLINE];

	while (gets(str) != NULL) {
		char *s = strchr(str, '(');
		s++;
		int n = strlen(s);
		for (int i = 0; i < n; ++i)
			if (s[i] == ',') s[i] = ' ';
		s[n - 1] = 0;
		int k;
		sscanf(s, "%d %s", &k, t);
		queryN++;
		queries.push_back(query(k, tag[string(t)]));
	}

	/* execute queries */
	pthread_t pt[THREADN];
	pthread_mutex_init(&mutex, NULL);

	for (int i = 0; i < THREADN; ++i) {
		pthread_create(&pt[i], NULL, doQuery4, 0);
	}

	for (int i = 0; i < THREADN; ++i) {
		pthread_join(pt[i], &status);
	}

	delete [] edge;
	delete [] p_hf_f;
	fclose(stdout);
	if(persons.size() > 6000000) {
		waitpid(query1_pid, NULL, 0);
	}
}

/**
 * Main query body
 */
void *doQuery4(void *pnt) {
	int k, t;
	while (1) {
		int total = 0;
		int current;
		/* read a query */
		pthread_mutex_lock(&mutex);
		if (curN == queryN) {
			pthread_mutex_unlock(&mutex);
			pthread_exit(NULL);
		}
		k = queries[curN].k;
		t = queries[curN].t;
		current = curN++;
		pthread_mutex_unlock(&mutex);

		/* run the query */
		runQuery(current, k, t);
	}
}

/**
 * Execute one query in the induced graph.
 */
void runQuery(int qid, int k, int t) {
	set <int> vforum;
	for (map <int, set<int> >::iterator iter = f_ht_t.begin(); iter != f_ht_t.end(); ++iter) {
		if ((*iter).second.find(t) != iter->second.end()) 
			vforum.insert((iter->first));
	}

	int n = persons.size();
	bool *vpl = new bool[n + 1]; /* valid person list */
	vector<int> valid_person_list;
	vector <int> *vedge = new vector <int>[n + 1]; // edge list
	memset(vpl, 0, sizeof(bool) * (n + 1));

	/* extract valid person set */
	int valid_person_size = 0;
	for (int person = 0; person < n; ++person) {
		bool key = 0;
		int size = p_hf_f[person].size();
		for (int i = 0; i< size; ++i)
			if (vforum.find(p_hf_f[person][i]) != vforum.end()) {
				key = 1;
				break;
			}
		if(!key) continue;
		vpl[person] = 1;
		valid_person_list.push_back(person);
		valid_person_size++;
	}

	/* extract graph here!*/
	for (int person = 0; person < n; ++person) {
		if (vpl[person]) {
			for (int t = 0; t < edge[person].size(); ++t) {
				int j = edge[person][t];
				if (vpl[j]) {
					vedge[person].push_back(j);
				}
			}
		}
	}
	delete[] vpl;
	vpl = NULL;

	/* find the top-k closensess centrality vertices */
	topkrank(qid, k, valid_person_list, vedge, valid_person_size, n);

	delete[] vedge;
}

/**
 * Process each connected component one by one.
 * 
 * In each connected component, using cyclic select-test method.
 */
void topkrank(int qid, int k, vector<int> valid_person_list, vector<int> *person_graph, int num_person, int total_num_person) {
	bool *visitedl = new bool[total_num_person + 1];
	vector<vertex_degree> cc_vertices;
	int cc_size = 0;
	queue <int> search_queue;
	memset(visitedl, 0, sizeof(bool) * (total_num_person + 1));
	set<vertex_centrality> answers;

	for (int idx=0; idx < valid_person_list.size(); ++idx) {
		int start = valid_person_list[idx];
		if (visitedl[start]) 
			continue;

		/* find one component using BFS. */
		search_queue.push(start);
		visitedl[start] = 1;
		cc_vertices.clear();
		cc_size = 0;
		while (!search_queue.empty()) {
			int cur = search_queue.front();
			search_queue.pop();
			cc_vertices.push_back(vertex_degree(person_graph[cur].size(), cur));
			++cc_size;
			for (int nbr = 0; nbr < person_graph[cur].size(); ++nbr) {
				int nv = person_graph[cur][nbr];
				if (visitedl[nv])
					continue;
				visitedl[nv] = 1;
				search_queue.push(nv);
			}
		}

		/* sort the vertices by the degree in ascending order.*/
		sort(cc_vertices.begin(), cc_vertices.end());

		/* execute K. Okamoto Algorithms */
		if(cc_size > 1)
			cyclic_calculation(qid, answers, k, cc_vertices, cc_size, person_graph, num_person, total_num_person);
		else if(cc_size == 1) {
			update_answer(answers, k, 0, persons[cc_vertices[0].id]);
		}
	}
	delete [] visitedl;

	/* output the answer. */
	char result[MAXLINE];
	char tmp[MAXTERM];
	sprintf(result, "%d", qid);
	for (set <vertex_centrality>::iterator iter = answers.begin(); iter != answers.end(); ++iter) {
		sprintf(tmp, " %ld", (*iter).id);
		strcat(result, tmp);
	}

	printf("%s\n", result);
}

/** 
 * Process procedure of delta_pfs.
 *
 */
void process(int i, int *dist, int sumDist, set <vertex_centrality> &answer, int *father, int k, vector <int> *vedge, int sumN, int count, bool *prunable, int &total, int level, short& is_changed) {
	level++;
	double old_centrality;

	if (sumN > 1) old_centrality = (double)(sumN - 1) * (sumN - 1) / sumDist / (count - 1);
	else old_centrality = 0;

	int n = persons.size();
	bool *v = new bool[n + 1];
	int *lambda = new int[n + 1];
	int beginT = clock();
	memmove(lambda, dist, sizeof(int) * (n + 1));
	int old_s = sumDist;
	for (int y = 0; y < vedge[i].size(); ++y) {
		int o = vedge[i][y];
		if (father[o] != i)
			continue;

		double centrality;
		if (sumN > 1) centrality = (double)(sumN - 1) * (sumN - 1) / (sumDist - sumN + 2) / (count - 1);
		else centrality = 0;
		if (answer.size() >= k && centrality < (answer.rbegin()->closeness - 1e-6) && prunable[o])
			continue;

		total++;

		sumDist += sumN - 2;
		dist[o] = dist[i] - 1;

		queue <int> list;
		memset(v, 0, sizeof(bool) * (n + 1));
		v[o] = 1;
		list.push(o);

		int cnt = 0;
		while (!list.empty()) {
			int cur = list.front();
			v[cur] = 1;
			list.pop();
			for (int t = 0; t < vedge[cur].size(); ++t) {
				int j = vedge[cur][t];
				cnt++;
				if (v[j])
					continue;
				v[j] = 1;
				if (dist[cur] + 1 >= dist[j])
					continue;
				list.push(j);
				sumDist -= dist[j] - (dist[cur] + 1);
				dist[j] = dist[cur] + 1;
			}
		}
		if (sumN > 1) centrality = (double)(sumN - 1) * (sumN - 1) / (sumDist) / (count - 1);
		else centrality = 0;
		answer.insert(vertex_centrality(centrality, persons[o]));
		if (answer.size() > k) {
			if(is_changed > 0 && answer.rbegin()->id != persons[o]) {
				is_changed *= -1;
			}
			answer.erase((*answer.rbegin()));
		}

		process(o, dist, sumDist, answer, father, k, vedge, sumN, count, prunable, total, level, is_changed);

		memmove(dist, lambda, sizeof(int) * (n + 1));

		sumDist = old_s;
	}
	delete[] lambda;
	delete[] v;
}

/**
 * Calculate exact clossness centraliry value here.
 *
 * Using the methon mentioned in "Efficient Top-k Closeness Centrality Search"
 * ICDE 2014
 */
void delta_pfs(set <vertex_centrality> &answer, int k, vector <vertex_degree> &vp, vector <int> *vedge, int total_num_person, int cc_size, int num_person, short& is_changed) {
	int n = total_num_person;

	bool *vpl = new bool[n + 1];
	bool *visitedl = new bool[n + 1];
	bool *v = new bool[n + 1];
	int *father = new int[n + 1];
	int *dist = new int[n + 1];
	bool *prunable = new bool[n + 1];
	queue <int> list;

	memset(vpl, 0, sizeof(bool) * (n + 1));
	memset(father, -1, sizeof(int) * (n + 1));
	memset(dist, 0, sizeof(int) * (n + 1));
	memset(prunable , 1, sizeof(bool) * (n + 1));

	for (int z = 0; z < vp.size(); ++z)
		vpl[vp[z].id] = 1;
	/*
	 * find a schedule and vertice number of each component
	 * */
	int total = 0;
	sort(vp.begin(), vp.end());	
	memset(visitedl, 0, sizeof(bool) * (n + 1));
	for (int z = 0; z < vp.size(); ++z) {
		int i = vp[z].id;
		if (visitedl[i]) continue;
		visitedl[i] = 1;
		int sumDist = 0, sumN = 0;		
		double centrality;

		memset(dist, 0, sizeof(int) * (n + 1));
		int maxd = 0;
		queue <int> schedule;
		schedule.push(i);

		/* calc dist map for vertex i */
		memset(v, 0, sizeof(bool) * (n + 1));
		v[i] = 1;
		list.push(i);
		while (!list.empty()) {
			int cur = list.front();
			v[cur] = 1;
			list.pop();
			sumDist += dist[cur];
			for (int t = 0; t < vedge[cur].size(); ++t) {
				int j = vedge[cur][t];
				if (v[j]) 
					continue;
				v[j] = 1;
				list.push(j);
				dist[j] = dist[cur] + 1;
			}
		}

		/* determine the schedular from vertex i */
		list.push(i);
		while (!list.empty()) {
			int cur = list.front();
			visitedl[cur] = 1;
			list.pop();
			for (int t = 0; t < vedge[cur].size(); ++t) {
				int j = vedge[cur][t];
				if (visitedl[j] || !vpl[j])
					continue;
				visitedl[j] = 1;
				list.push(j);
				father[j] = cur;
				prunable[cur] = 0;
			}
		}

		if (cc_size > 1) centrality = (double)(cc_size - 1) * (cc_size - 1) / (sumDist) / (num_person - 1);
		else centrality = 0;
		answer.insert(vertex_centrality(centrality, persons[i]));
		if (answer.size() > k) {
			if(is_changed > 0 && answer.rbegin()->id != persons[i]) {
				is_changed *= -1;
			}
			answer.erase((*answer.rbegin()));
		}

		/* from now, begin the delta-pfs */
		process(i, dist, sumDist, answer, father, k, vedge, cc_size, num_person, prunable, total, 0, is_changed);
	}
	delete [] vpl;
	delete [] visitedl;
	delete [] v;
	delete [] father;
	delete [] dist;
	delete [] prunable;
}



/**
 * cyclic_calculation
 * 
 * Select some candidates by approximate method, calculate the exact centraliry value and update the top-k answers.
 * If the top-k answer set changed, do it again!
 */
void cyclic_calculation(int qid, set<vertex_centrality>& answers, int k, vector<vertex_degree>& cc_vertices, int cc_size, vector<int> *person_graph, int num_person, int total_num_person) {
	int *exact2hop = new int[total_num_person + 1];

	int num_samples = 0;
	int topke = 200 * k; // choose top 100*k of estimated closeness
	if (total_num_person > 600000 && k < 5) topke = 1000;
	int topd = 0; // choose top 100 biggest of degree
	double alpha = 1; /* tunning later. */

	vector<vertex_avgdist> estimated_avgdist;
	vector<int> candidates;
	int* estimated_distances = new int[total_num_person+1];
	bool* is_in_candidates = new bool[total_num_person+1];
	int estimated_radius = -1; /**/
	memset(estimated_distances, 0, sizeof(int)*(total_num_person+1));
	memset(is_in_candidates, 0, sizeof(bool)*(total_num_person+1));

	/* using heuristic rules to select seeds and estimate the centrality! */
	for(int i = 0; i < num_samples; ++i) {
		int start = cc_vertices[i].id;
		int tmp_radius = -1;
		bfs(start, estimated_distances, person_graph, total_num_person, tmp_radius);
		if(tmp_radius < estimated_radius || estimated_radius == -1)
			estimated_radius = tmp_radius;
	}

	int n = total_num_person;
	long **sketch[2];	
	int *hashtable = new int[n + 1];
	int *total_estimated = new int[n + 1];
	memset(total_estimated, 0, sizeof(int) * (n + 1));
	sketch[0] = new long*[n + 1];
	sketch[1] = new long*[n + 1];
	int mask = NMAP - 1;
	int offset = popcnt(NMAP - 1);
	int bigmask = (1 << 16) - 1;
	int *hopin2 = new int[n + 1];
	memset(hopin2, -1, sizeof(int) * (n + 1));
	memset(exact2hop, 0, sizeof(int) * (n + 1));
	for (int iter = 0; iter < cc_size; ++iter) {
		int i = cc_vertices[iter].id;
		for (int j = 0; j < person_graph[i].size(); ++j) {
			int cur = person_graph[i][j];
			if (hopin2[cur] < iter) {
				hopin2[cur] = iter;
				exact2hop[i]++;
			}
			int length = person_graph[cur].size();
			for (int z = 0; z < length; ++z) {
				int u = person_graph[cur][z];
				if (hopin2[u] < iter) {
					hopin2[u] = iter;
					exact2hop[i]++;
				}

			}
		}
		estimated_distances[i] -= person_graph[i].size();
		estimated_distances[i] -= exact2hop[i];
	}

	delete [] hopin2;
	delete [] hashtable;
	for (int i = 0; i < 2; ++i) {
		delete [] sketch[i];
	}


	for(int i = 0; i < cc_size; ++i) {
		int vid = cc_vertices[i].id;
		double avgdist = estimated_distances[vid];
		estimated_avgdist.push_back(vertex_avgdist(avgdist, vid));
	}
	sort(estimated_avgdist.begin(), estimated_avgdist.end());

	int bound = 1 << 28;
	if (cc_size >= k) {
		bound = estimated_avgdist[k - 1].avgdist + (double)total_estimated[estimated_avgdist[0].id] * 1.65 * 0.068 * 0.4;
	}

	vector <vertex_degree> vp;
	/* calculate the candidate E */
	for(int i = 0; i < topke && i < estimated_avgdist.size(); ++i) {
		candidates.push_back(estimated_avgdist[i].id);
		is_in_candidates[estimated_avgdist[i].id] = 1;
		vp.push_back(vertex_degree(exact2hop[estimated_avgdist[i].id], estimated_avgdist[i].id));
	}
	for(int i = cc_vertices.size() - 1; i + topd >= cc_vertices.size() && i >=0; --i) {
		if(is_in_candidates[cc_vertices[i].id] == 0) {
			candidates.push_back(cc_vertices[i].id);
			vp.push_back(vertex_degree(cc_vertices[i].degree, cc_vertices[i].id));
		}
	}

	delete [] is_in_candidates;
	estimated_distances = NULL;
	is_in_candidates = NULL;

	/* find the top-k in this connected component and merged into answer. */
	short is_changed = 0;
	int inc = 70*k;
	delta_pfs(answers, k, vp, person_graph, total_num_person, cc_size, num_person, is_changed);
	do{
		is_changed *= -1;
		is_changed++;
		vp.clear();
		for(int i = topke; i < topke+inc && i < estimated_avgdist.size(); ++i) {
			vp.push_back(vertex_degree(exact2hop[estimated_avgdist[i].id], estimated_avgdist[i].id));
		}
		topke += inc;
		if(vp.size() > 0) {
			delta_pfs(answers, k, vp, person_graph, total_num_person, cc_size, num_person, is_changed);
		}
	}while(is_changed < 0 && vp.size() > 0);
	delete [] exact2hop;
	estimated_avgdist.clear();
}

void update_answer(set<vertex_centrality>& answer, int k, double centrality, long vid) {
	answer.insert(vertex_centrality(centrality, vid));
	if (answer.size() > k) 
		answer.erase((*answer.rbegin()));
}

void bfs(int start, int* estimated_distances, vector <int> * person_graph, int total_num_person, int& estimated_radius) {
	int *dist = new int[total_num_person + 1];
	queue <int> search_queue;
	memset(dist, -1, sizeof(int) * (total_num_person + 1));

	/* compute shortest path using BFS. */
	search_queue.push(start);
	dist[start] = 0;

	while(!search_queue.empty()) {
		int cur = search_queue.front();
		search_queue.pop();
		estimated_distances[cur] += dist[cur];
		if(dist[cur] > estimated_radius)
			estimated_radius = dist[cur];

		for(int nbr = 0; nbr < person_graph[cur].size(); ++nbr) {
			int nv = person_graph[cur][nbr];
			if (dist[nv] != -1)
				continue;
			dist[nv] = dist[cur] + 1;
			search_queue.push(nv);
		}
	}
	delete [] dist;
}

