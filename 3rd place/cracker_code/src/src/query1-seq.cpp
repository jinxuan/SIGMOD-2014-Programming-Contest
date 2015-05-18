/**
 * Sequential version of Query Type 1.
 *
 * Detailed explaintions can be found in query1.cpp.
 */

#include <sys/syscall.h>
#include <unistd.h>
#include <sched.h>

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>
#include <vector>
#include <io_mmap.hpp>
#include <google/dense_hash_map>
#include <sparsehash/sparse_hash_map>
#include <meminfo.hpp>

#define gettid() syscall(__NR_gettid) 
using namespace std;
using std::tr1::hash;

#define CONCURRENT_UPDATE
#define THREADN 8
#define PQUERY_THRESHOLD 60000

char dir[100] = "";
int n_max = 10000; /* record the maximal number of persons. */
int m_max = 10000;
long c_max = 1000000;
int oo = 1000000000;

int* cnt; /* record the comment count. */
vector<int> src;
vector<int> dst;

vector<int> tmp_src;
vector<int> tmp_dst;

/* store the comment_hasCreator_person */
//vector<long> comments;
vector<int> creators;
int n, m;

/* remap the person id */
struct eqlong {
	bool operator()(const long a, const long b) const {
		return a == b;
	}
};
google::dense_hash_map<long, int, hash<int>, eqlong> p2i;   //person id --> new id
google::sparse_hash_map<long, int, hash<long>, eqlong> c2p;

bool is_use_sparse = false;

/* global reader */
mmap_reader mmap_file_reader;

/* adjancent list */
typedef struct Nbs {
	int n, head, tail;
	int find(int p) {
		int l = head, r = tail-1;
		while (l <= r) {
			int m = (l+r)/2;
			if (dst[m] == p) return cnt[m];
			if (p < dst[m]) r = m-1;
			else l = m+1;
		}
		return oo;
	}
	int up(int p) {
		int l = head, r = tail-1;
		while (l <= r) {
			int m = (l+r)/2;
			if (dst[m] == p) {
				__sync_add_and_fetch(cnt+m, 1);
				return 1;
			}
			if (p < dst[m]) r = m-1;
			else l = m+1;
		}
		return 0;
	}
	int change(int p, int x) {
		int l = head, r = tail-1;
		while (l <= r) {
			int m = (l+r)/2;
			if (dst[m] == p) {
				if (cnt[m]>x) cnt[m] = x;
				return 1;
			}
			if (p < dst[m]) r = m-1;
			else l = m+1;
		}
		return 0;
	}
	void qsort(int l, int r) {
		if (l>=r) return;
		int i = l, j = r, x = cnt[(l+r)/2];
		while (i<=j) {
			while (i<=j && cnt[i]>x) ++i;
			while (i<=j && cnt[j]<x) --j;
			if (i <= j) {
				int tmp = dst[i]; dst[i] = dst[j]; dst[j] = tmp;
				tmp = cnt[i]; cnt[i] = cnt[j]; cnt[j] = tmp;
				++i;--j;
			}
		} 
		if (l<j) qsort(l, j);
		if (i<r) qsort(i, r);
	}
	void sort() {
		qsort(head, tail-1);
	}
}Nbs;

Nbs* f;

void bsort_src(int n, int m) {
	for (int i = 0; i <= n; ++i) cnt[i] = 0;
	for (int i = 0; i < m; ++i) cnt[tmp_src[i]]++;
	for (int i = 1; i <= n; ++i) cnt[i] += cnt[i-1];
	for (int i = m-1; i >= 0; --i) {
		int k = --cnt[tmp_src[i]];
		src[k] = tmp_src[i]; dst[k] = tmp_dst[i];
	}
}

void bsort_dst(int n, int m) {
	for (int i = 0; i <= n; ++i) cnt[i] = 0;
	for (int i = 0; i < m; ++i) cnt[dst[i]]++;
	for (int i = 1; i <= n; ++i) cnt[i] += cnt[i-1];
	tmp_src.reserve(dst.size());
	tmp_dst.reserve(dst.size());
	for (int i = m-1; i >= 0; --i) {
		int k = --cnt[dst[i]];
		tmp_src[k] = src[i]; tmp_dst[k] = dst[i];
	}
}

void bsort(int n, int m) {
	bsort_dst(n, m);
	bsort_src(n, m);
}

long max_id; /* record the max id for estimating the number of lines in person_knows_person */
void getperson(char *dbpath) {
	strcpy(dir, dbpath);

	p2i.set_empty_key(-1);
	p2i.resize(n_max);
	n = 0;
	max_id = 0;

	bool success = mmap_file_reader.open_mmap(strcat(dir,"/person.csv"));
	if(success) {
		mmap_file_reader.skipfirstline();

		long id;
		while (mmap_file_reader.hasNext()) {
			id = mmap_file_reader.getLong();
			if(max_id < id) max_id = id;
			p2i[id] = n++;
			mmap_file_reader.skipline();
		}
		mmap_file_reader.close_mmap();
	}
	else{
		freopen(dir, "r", stdin);
		char str[1024];
		char id[100];
		char* s;
		long x;
		gets(str);
		while (gets(str) != NULL) {
			int i = 0;
			s = str;
			while (*s != '|') id[i++] = *(s++);
			s++;
			id[i] = 0;
			x = strtol(id, NULL, 10);
			if(max_id < x) max_id = x;
			p2i[x] = n++;
		}
	}
}

int getdigitalcount(long num) {
	int cnt = 0;
	while(num > 0) {
		++cnt;
		num /= 10;
	}
	return cnt > 1 ? cnt - 1 : 0;
}

void init1(char* dbpath) {
	getperson(dbpath);

	strcpy(dir, dbpath);
	bool success = mmap_file_reader.open_mmap(strcat(dir, "/person_knows_person.csv"));

	m = 0;
	long a,b;
	long fs = mmap_file_reader.filesize();
	int unit = getdigitalcount(n - 1) * 2 + 2 + getdigitalcount((max_id+n-1)/n);
	m_max = fs / unit;
	src.reserve(m_max);
	dst.reserve(m_max);
	if(success) {
		mmap_file_reader.skipfirstline();
		while (mmap_file_reader.hasNext()) {
			a = mmap_file_reader.getLong();
			b = mmap_file_reader.getLong();
			src.push_back(p2i[a]);
			dst.push_back(p2i[b]);
			m++;
		}
		mmap_file_reader.close_mmap();
	}
	else{
		char str[1024];
		char id[100], id2[100];
		char* s;
		freopen(dir, "r", stdin);
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
			src.push_back(p2i[strtol(id, NULL, 10)]);
			dst.push_back(p2i[strtol(id2, NULL, 10)]);
			m++;
		}
	}

	f = new Nbs[n_max];
	cnt = new int[n_max > m ? n_max : m]; /* use to record the comment count later!! */
	bsort(n, m);

	int head = 0, tail = 0;
	for (int i = 0; i < n; ++i) {
		while (tail < m && src[tail] == i) {
			++tail;
		}
		f[i].head = head; f[i].tail = tail; f[i].n = tail-head;
		head = tail;
	}

	/* free memory */
	tmp_src.clear();
	tmp_dst.clear();
	src.clear();
}

int delta;

void init2(char* dbpath) {
	strcpy(dir, dbpath);
	bool success = mmap_file_reader.open_mmap(strcat(dir,"/comment_hasCreator_person.csv"));

	/* estimated the c_max */
	long fs = mmap_file_reader.filesize();
	int unit;
	c_max = n * 200;
	if(n < 2000) unit = 9;
	else if(n < 20000) unit = 11;
	else if(n < 200000) unit = 13;
	else unit = 16;
	fs = fs / unit;
	c_max = fs > c_max ? fs : c_max;

	if(success) {

		long offset = 0, readbytes = 0;

		readbytes = mmap_file_reader.skipfirstline();
		creators.reserve(c_max);
		int cur_max = c_max;
		long a, b; 

		long limit = 5L * 1024 * 1024 * 1024; //5G

		/* reader first record*/
		readbytes += mmap_file_reader.getLong(a);
		readbytes += mmap_file_reader.getLong(b);
		if((a /delta) > cur_max) {
			cur_max = a / delta;
			creators.resize(cur_max + 1);
		}
		creators[a / delta]=p2i[b];

		while (mmap_file_reader.hasNext()) {
			if(readbytes > limit) {
				offset += readbytes;
				readbytes = 0;
				mmap_file_reader.remap(dir, offset);
			}
			readbytes += mmap_file_reader.getLong(a);
			readbytes += mmap_file_reader.getLong(b);
			if((a /delta) > cur_max) {
				cur_max = a / delta;
				creators.resize(cur_max + 1);
			}
			creators[a/delta] = p2i[b];
		}
		mmap_file_reader.close_mmap();
	}
	else{
		freopen(dir, "r", stdin);
		char str[1024];

		int cur_max = c_max;

		creators.reserve(c_max);

		gets(str);
		long pre = -10000000, a, b;
		scanf("%ld|%ld\n", &a, &b);

		if((a /delta) > cur_max) {
			cur_max = a / delta;
			creators.resize(cur_max + 1);
		}
		creators[a / delta]=p2i[b];


		while (scanf("%ld|%ld\n", &a, &b) != EOF) {
			if((a /delta) > cur_max) {
				cur_max = a / delta;
				creators.resize(cur_max + 1);
			}
			creators[a/delta] = p2i[b];
		}
	}
}

long file_block_size;
void* concurrent_aggregate(void* start_) {
	long start = *(long*)start_;
	if(n < PQUERY_THRESHOLD) {
		cpu_set_t *cpusetp;
		cpusetp = CPU_ALLOC(4);
		CPU_SET(0, cpusetp);
		CPU_SET(2, cpusetp);
		CPU_SET(4, cpusetp);
		CPU_SET(6, cpusetp);
		if(sched_setaffinity(gettid(), sizeof(cpusetp), cpusetp) == -1) {}
	}

	long has_read = mmap_file_reader.skiplinefrom(start);
	long p1, p2;
	while(has_read <= file_block_size && start < mmap_file_reader.map_file_size()) {
		/* read a line */
		has_read += mmap_file_reader.getLong(start, p1);
		has_read += mmap_file_reader.getLong(start, p2);
		f[creators[p1/delta]].up(creators[p2/delta]);
	}
}

void init3(char* dbpath) {
	strcpy(dir, dbpath);
	strcat(dir,"/comment_replyOf_comment.csv");

	for (int i = 0; i < m; ++i)
		cnt[i] = 0;//0 -> cnt

	if(mmap_file_reader.split_file(dir)) {
		int blk_num = mmap_file_reader.get_fileblock_num();
		for(int k = 0; k < blk_num; ++k) {
			mmap_file_reader.open_mmap(dir, k);

			pthread_t pt[THREADN];
			long offsets[THREADN];

			file_block_size = (mmap_file_reader.fb[k].length -mmap_file_reader.fb[k].shift+ THREADN - 1) / THREADN;
			for(int i = 0; i < THREADN; ++i) {
				offsets[i] = file_block_size * i + mmap_file_reader.fb[k].shift;
			}

			for (int i = 0; i < THREADN; ++i) {
				pthread_create(&pt[i], NULL, concurrent_aggregate, (void*)(offsets+i));
			}

			for (int i = 0; i < THREADN; ++i) {
				pthread_join(pt[i], NULL);
			}
			mmap_file_reader.close_mmap();
		}
	}
	else{
		mmap_file_reader.open_mmap(dir);
		mmap_file_reader.skipfirstline();
		/* parallel execute bfs */
		pthread_t pt[THREADN];
		long offsets[THREADN];

		file_block_size = (mmap_file_reader.filesize() + THREADN - 1) / THREADN;
		for(int i = 0; i < THREADN; ++i) {
			offsets[i] = file_block_size * i;
		}

		for (int i = 0; i < THREADN; ++i) {
			pthread_create(&pt[i], NULL, concurrent_aggregate, (void*)(offsets+i));
		}

		for (int i = 0; i < THREADN; ++i) {
			pthread_join(pt[i], NULL);
		}
		mmap_file_reader.close_mmap();
	}

	creators.clear();

	for (int i = 0; i < n; ++i)
		for (int l = f[i].head; l < f[i].tail; ++l)
			f[dst[l]].change(i, cnt[l]);
	for (int i = 0; i < n; ++i)
		f[i].sort();
}

int* dist;
int* dq;
int* v;
bool* tag; //check p2
int vs = 0;

int query1(long p1_, long p2_, int x) {
	if (p1_ == p2_) return 0;

	int p1 = p2i[p1_];
	int p2 = p2i[p2_];

	++vs;
	dist[p1] = 0;
	dist[p2] = 0;
	v[p1] = vs;
	v[p2] = vs;
	tag[p1] = false;
	tag[p2] = true;
	dq[0] = p1;
	dq[1] = p2;

	int l = 0, r = 2;
	while (l < r) {
		int p = dq[l];
		for (int i = f[p].head; i < f[p].tail; ++i) {
			int q = dst[i];
			if (cnt[i]>x) {
				if ((v[q]==vs) && (tag[p]^tag[q])) {
					return dist[p]+dist[q]+1;
				} else {
					if (v[q]!=vs) {
						dist[q] = dist[p] + 1;
						tag[q] = tag[p];
						v[q] = vs;
						dq[r++] = q;
					}
				}
			} else break;
		}
		++l;
	}
	return -1;
}

int proc(char* query_path, char* out_path) {
	freopen(query_path, "r", stdin);
	freopen(out_path, "w", stdout);
	long p1, p2;
	int x;
	int qid = 0;

	dist = new int[n_max];
	dq = new int[n_max];
	v = new int[n_max];
	tag = new bool[n_max];

	memset(v,0,sizeof(int)*n_max);

	while(scanf("query1(%ld, %ld, %d)\n", &p1, &p2, &x)!=EOF) {
		printf("%d %d\n", qid, query1(p1, p2, x));
		++qid;
	}
	delete[] dist;
	delete[] dq;
	delete[] v;
	delete[] tag;
	delete[] f;
	delete[] cnt;
	return qid;
}

/* parallel query processing */
class query {
	public:
		int p1, p2, x;
		query(int p1_, int p2_, int x_):p1(p1_), p2(p2_), x(x_) {
		}
};
vector<query>* queries;

void *doQuery1(void* pnt) {
	int tid = *(int*)pnt;
	int p1, p2, x;
	int* dist = new int[n_max];
	int* dq = new int[n_max];
	int* v = new int[n_max];
	bool* tag = new bool[n_max];
	memset(v, 0, sizeof(int)*n_max);
	int vs = 0;

	int size = queries[tid].size();
	for(int i = 0; i < size; i++) {
		int current = i * THREADN + tid;
		p1 = queries[tid][i].p1;
		p2 = queries[tid][i].p2;
		x = queries[tid][i].x;

		/* do the query! */
		int ans = -1;
		if (p1 == p2) { ans  = 0; }
		else{
			++vs;
			dist[p1] = 0;
			dist[p2] = 0;
			v[p1] = vs;
			v[p2] = vs;
			tag[p1] = false;
			tag[p2] = true;
			dq[0] = p1;
			dq[1] = p2;

			int l = 0, r = 2;
			while (l < r && ans == -1) {
				int p = dq[l];
				for (int i = f[p].head; i < f[p].tail; ++i) {
					int q = dst[i];
					if (cnt[i]>x) {
						if ((v[q]==vs) && (tag[p]^tag[q])) {
							ans = dist[p]+dist[q]+1;
							break;
						} else {
							if (v[q]!=vs) {
								dist[q] = dist[p] + 1;
								tag[q] = tag[p];
								v[q] = vs;
								dq[r++] = q;
							}
						}
					} else break;
				}
				++l;
			}
		}
		printf("%d %d\n", current, ans);
	}

	delete[] dist;
	delete[] dq;
	delete[] v;
	delete[] tag;
}

int pproc(char* query_path, char* out_path) {
	freopen(query_path, "r", stdin);
	freopen(out_path, "w", stdout);
	long p1, p2;
	int x;

	queries = new vector<query>[THREADN];
	int qid = 0;

	while(scanf("query1(%ld, %ld, %d)\n", &p1, &p2, &x)!=EOF) {
		queries[qid % THREADN].push_back(query(p2i[p1], p2i[p2], x));
		++qid;
	}

	pthread_t pt[THREADN];

	for (int i = 0; i < THREADN; ++i) {
		pthread_create(&pt[i], NULL, doQuery1, new int(i));
	}

	for (int i = 0; i < THREADN; ++i) {
		pthread_join(pt[i], NULL);
	}
	delete[] f;
	delete[] cnt;
	return 0;
}


int main(int argc, char** argv) {
	time_t start_time, end_time;
	start_time = time(NULL);

	n_max = atoi(argv[4]);
	delta = atoi(argv[5]);

	init1(argv[1]);

	init2(argv[1]);

	init3(argv[1]);

	if(n < PQUERY_THRESHOLD)
		proc(argv[2], argv[3]);
	else
		pproc(argv[2], argv[3]);
}
