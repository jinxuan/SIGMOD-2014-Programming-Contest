#include <cstdio>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <cstring>
#include <string>
#include <set>
#include <vector>
#include <functional>

using namespace std;

#include <io_mmap.hpp>
#include <google/dense_hash_map>
#include <meminfo.hpp>

using google::dense_hash_map;
using std::tr1::hash;

#define MAXLINE 500
#define MAXTERM 1024

/** initial IO related */
char dir[MAXTERM];

struct Cmpint {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

struct Cmp2int {
	bool operator()(const int &a, const int &b) {
		return a > b;
	}
};

struct Cmp3{
	bool operator() (const pair<int, string> &a, const pair<int, string> &b) {
		if (a.first != b.first) return a.first > b.first;
		else return a.second < b.second;
	} 
};

struct eqlong {
	bool operator()(const long a, const long b) const {
		return a == b;
	}
};
struct eqint {
	bool operator()(const int a, const int b) const {
		return a == b;
	}
};

google::dense_hash_map<int, set<int>, hash<int>, eqint> person_has_tag;   /** person id --> tag id */
map <int, string>  tag;  /** tag(id,name), maintain a one-one map between id and name */
map <pair<int, int>, int, Cmpint> query2;  /** record queries, sort by date and maintain the original order */
multimap <int, int, Cmp2int> personinfo; /** person.birth <-> person.id (one-one map & sort by date) */
vector<vector<string> > ans;  /** record the answer. ans[i][j] represents the (j + 1)th intereset tag in the (i+1)th query */
google::dense_hash_map<int, int, hash<int>, eqint>ans1; /** tag range number correspondent to every query. pair<int,int>  <-> (tag.id, tag.range) */
vector<int>* edge;  /** adjacent list of person knows person */
google::dense_hash_map<long, int, hash<int>, eqlong> p2i;    /** person id --> new id */

int n_max;

int pn = 0, m = 0, person_has_tag_size =0; // people number & query number

/* declare for disjoint set */
void disjoint_set_init(int person_hasinterest_tag_size);
void disjoint_set_clear();
void add_disjoint_set_node(int tag, int pid);
int find(int tag, int pid);
int merge(int cc1, int cc2);

int convert_date_to_int(string date) {
	int ans = 0;
	for(string::iterator it = date.begin(); it != date.end(); ++it){
		if(*it == '-') continue;
		ans = ans * 10 + (*it) - 48;
	}
	return ans;
}

void getperson(char *dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;
	if(mmap_file_reader.open_mmap(strcat(dir,"/person.csv")) == false)
		return;
	mmap_file_reader.skipfirstline();

	long id;
	char birth[MAXTERM];

	p2i.set_empty_key(-1);
	p2i.resize(n_max);

	pn = 0;
	while (mmap_file_reader.hasNext()){
		id = mmap_file_reader.getLong();
		mmap_file_reader.skiprecord(); //first name
		mmap_file_reader.skiprecord(); //last name
		mmap_file_reader.skiprecord(); //gender
		mmap_file_reader.getString(birth);

		personinfo.insert(make_pair(convert_date_to_int(string(birth)), pn));
		p2i[id] = pn++;

		mmap_file_reader.skipline();
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
		tag[id] = string(name);
		mmap_file_reader.skipline();
	}
	mmap_file_reader.close_mmap();
}

void getperson_hasInterest_tag(char *dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;
	if(mmap_file_reader.open_mmap(strcat(dir,"/person_hasInterest_tag.csv")) == false)
		return;
	mmap_file_reader.skipfirstline();

	long id;
	int id2;
	person_has_tag.set_empty_key(-1);
	person_has_tag.resize(pn + 1);
	while(mmap_file_reader.hasNext()) {
		id = mmap_file_reader.getLong();
		id2 = mmap_file_reader.getInt();
		person_has_tag[p2i[id]].insert(id2);
		person_has_tag_size++;
	}
	mmap_file_reader.close_mmap();
}

/**
 *  Read the graph.
 *  Must be called after calling getperson method!!!
 */
void getperson_knows_person(char *dbpath) {
	strcpy(dir, dbpath);
	mmap_reader mmap_file_reader;
	if(mmap_file_reader.open_mmap(strcat(dir,"/person_knows_person.csv")) == false)
		return;
	mmap_file_reader.skipfirstline();

	int id, id2;
	edge = new vector<int>[pn+1];
	int n = pn; /* here need the getperson method has been called. */
	for (int i = 0; i < n; ++i)
		edge[i].reserve(100);
	while(mmap_file_reader.hasNext()) {
		id = mmap_file_reader.getLong();
		id2 = mmap_file_reader.getLong();
		edge[p2i[id]].push_back(p2i[id2]);
	}
	mmap_file_reader.close_mmap();
}

/**
 * Function method
 *
 * @t Interest tag.
 * @pid Person ID.
 * @return Whether person with pid has interest tag t.
 */
bool check(int t, int pid) {
	if (person_has_tag[pid].find(t) != person_has_tag[pid].end())
		return true;
	return false;
}

bool ok(pair<int, string> temp, google::dense_hash_map<int, int, hash<int>, eqint>::iterator b) {
	if (temp.first < b->second) return true;
	if (temp.first == b->second && temp.second > tag[b->first]) return true;
	return false;
}

void getans(int qid, int k) {
	map<pair<int, string>, int, Cmp3> update;
	int n = ans1.size();
	google::dense_hash_map<int, int, hash<int>, eqint>::iterator partial_ans_it = ans1.begin();
	for (int j = 0; j < k; ++j) {
		pair<int, string> temp = make_pair(partial_ans_it->second, tag[partial_ans_it->first]);
		update[temp] = j;
		++partial_ans_it;
	}  

	for (int j = k; j < n && partial_ans_it != ans1.end(); ++j) {
		pair<int, string> uu = update.rbegin() -> first;
		if (ok(uu, partial_ans_it)) {
			map<pair<int, string>, int, Cmp3>::reverse_iterator ur_Iter = update.rbegin();         
			update.erase(ur_Iter -> first);
			pair<int, string> temp = make_pair(partial_ans_it->second, tag[partial_ans_it->first]);
			update[temp] = j;                                                                                                
		}
		++partial_ans_it;
	}
	ans[qid].resize(k);
	map<pair<int, string>, int, Cmp3>::iterator u_Iter;
	int r = 0;
	for (u_Iter = update.begin(); u_Iter != update.end(); u_Iter++) {
		ans[qid][r] = u_Iter -> first.second;
		r++;
	}
}

/**
 * Here to execute the calculation of Query Type 2.
 *
 * Use find-union algorithm. Add person incrementally by date. Collect answer during the process. 
 */
void doQuery2() {
	map<pair<int, int>, int, Cmpint>::iterator s_Iter; /** query iterator */
	map <int, int, Cmp2int>::iterator p_Iter; /** person itrator */
	int date_now; /** date of current query */
	map<int, vector<int> > size;
	bool* now;

	/* for every tag, update by date */
	s_Iter = query2.begin();
	p_Iter = personinfo.begin();
	now = new bool[pn+1]; /** record the valid person */
	memset(now, 0, sizeof(bool)*(pn+1));

	disjoint_set_init(person_has_tag_size);

	ans1.set_empty_key(-1024);
	ans1.resize(tag.size() + 1);
	for(map<int, string>::iterator t_Iter = tag.begin(); t_Iter != tag.end(); ++t_Iter){
		ans1[t_Iter->first] = 0;
	}
	while (s_Iter != query2.end() ) {
		int tag_num = 0;
		int nb_num = 0;

		date_now = s_Iter -> first.first;
		int max = 0;

		for( ; p_Iter != personinfo.end(); ++p_Iter){
			if(p_Iter->first >= date_now){
				// process_person++;
				now[p_Iter->second] = 1;

				/* update tag here */
				/* Disjoint set and find-union algoirthm. */
				set<int> person_tags = person_has_tag[p_Iter->second];

				for (set<int>::iterator tag_it = person_tags.begin(); tag_it != person_tags.end(); ++tag_it ){
					int id_ = p_Iter->second;
					int nb = edge[id_].size();

					add_disjoint_set_node(*tag_it, id_);
					max = 1;

					for(int i = 0 ; i < nb; ++i){
						if(now[edge[id_][i]] && check(*tag_it, edge[id_][i])){
							nb_num++;
							int x = edge[id_][i];
							int l = find(*tag_it, x);
							int j = find(*tag_it, id_);
							int newsize = -1;
							if(l != j){
								newsize = merge(l,j);
								if (newsize > max) 
									max = newsize;
							}
						}
					}

					google::dense_hash_map<int, int, hash<int>, eqint>::iterator tmpit = ans1.find(*tag_it);
					if(tmpit == ans1.end()){
						ans1[*tag_it] = max;
					}
					else if(tmpit->second < max)
						tmpit->second = max;
				}
			}
			else break;
		}

		getans(s_Iter->second, s_Iter->first.second);

		s_Iter++;
	}
	delete [] now;
	disjoint_set_clear();
}

int main(int argc, char** argv) {
	/* initialization. */ 
	char* dbpath = argv[1];
	n_max = atoi(argv[4]);

	getperson(dbpath);
	gettag(dbpath);
	getperson_hasInterest_tag(dbpath);
	getperson_knows_person(dbpath);

	/* Record all queries. Sort queries by date using map's property. */ 
	freopen(argv[2], "r", stdin);
	freopen(argv[3], "w", stdout);
	m = 0; /* m is the query number */
	int yy, mm, dd, k;
	while(scanf("query2(%d, %d-%d-%d)\n", &k, &yy, &mm, &dd) != EOF){
		pair<int, int> temp = make_pair(yy*10000+mm*100+dd, k);
		query2[temp] = m;
		m++; 
	}
	ans.resize(m);

	/* running query2 */
	doQuery2();

	fclose(stdin);

	/* output answer */ 
	for (int i = 0; i < m; i++) {
		printf("%d", i);
		for (int j = 0; j < ans[i].size(); j++)
			printf(" %s", ans[i][j].c_str());
		printf("\n");
	}
}

/******************************************************* Here are the variables and functions for find-union algorithm. **********************************************************************/

int* father;
int* size;
int used_space;

/* multi tag-related disjoint sets. */
map <int, map<int, int> > multi_disjoint_sets; /** each tag maintians a disjoint set. */

void disjoint_set_init(int person_hasinterest_tag_size){
	father = new int[person_hasinterest_tag_size + 10];
	size = new int[person_hasinterest_tag_size + 10];
	used_space = 0;
}

void disjoint_set_clear(){
	delete [] father;
	delete [] size;
	used_space = 0;
}

void add_disjoint_set_node(int tag, int pid){
	int id = used_space++;
	father[id] = id;
	size[id] = 1;
	multi_disjoint_sets[tag][pid] = id;
}

int internal_find(int cc){
	if(father[cc] != cc) 
		father[cc] = internal_find(father[cc]);
	return father[cc];
}

int find(int tag, int pid) {
	int cc = multi_disjoint_sets[tag][pid];
	return internal_find(cc);
}

/**
 * Merge two disjoint sets.
 *
 * @cc1 The first set.
 * @cc2 The second set.
 * @return Size of the merged set.
 */
int merge(int cc1, int cc2){
	/* merge small tree into large tree. */
	if (size[cc1] > size[cc2]) {
		father[cc2] = cc1;
		size[cc1] += size[cc2];
		return size[cc1];
	}
	else{
		father[cc1] = cc2;
		size[cc2] += size[cc1];
		return size[cc2];
	} 
}
