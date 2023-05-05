#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <future>
#include <thread>
#include "gif-h/gif.h"
using namespace std;

//ideally we want the dimensions of the field to be divisible by 60 or even 120
const int N = 1200; 
const int M = 1800;

typedef unsigned int uint;

uint matrix[N][M];
bool updated[N][M];

mt19937 random_engine;

const uint CELLTYPE_MASK = 3;
const uint REPRODUCTION_COUNTER_MASK = 4+8+16+32;
const uint ENERGY_COUNTER_MASK = 128+256;
const uint SHARK_LIFETIME = 4;
const uint SHARK_REPRODUCTION_TIME = 5;
const uint HERING_REPRODUCTION_TIME = 2;
enum CellType { empty=0, hering=1, shark=2 };

uint get_celltype(uint x) { return x&CELLTYPE_MASK; }
uint get_reproduction_counter(uint x) { return x&REPRODUCTION_COUNTER_MASK; }
uint get_energy_counter(uint x) { return x&ENERGY_COUNTER_MASK; }
uint set_celltype(uint x, CellType type) { return (x&~CELLTYPE_MASK) | type; }
uint increment_counter(uint x, uint mask, uint lifetime) {
	const uint mask_first_bit = (mask&-mask);
	return (x&~mask) | ((((x&mask) / mask_first_bit + 1) % lifetime) * mask_first_bit);
}

uint clear_energy_counter(uint x) { return x&~ENERGY_COUNTER_MASK; }
uint clear_reproduction_counter(uint x) { return x&~REPRODUCTION_COUNTER_MASK; }

void visualize_console_frame() {
	for (int i=0; i<N; ++i) {
		for (int j=0; j<M; ++j) {
			auto m = matrix[i][j];
			switch(get_celltype(m)) {
				case empty:  cout << "0"; break;
				case hering: cout << "*"; break;
				case shark:  cout << "#"; break;
			}
		}
		cout << endl;
	}
}

void visualize_frame(vector<uint8_t> &frame) {
	for (int p=0; p<frame.size(); p += 4) {
		int i = p/4;
		switch (get_celltype(matrix[i/M][i%M])) {
			case shark:  frame[p] = 220; frame[p+1] = 220; frame[p+2] = 220; frame[p+3] = 255; break;
			case hering: frame[p] = 220; frame[p+1] = 0;   frame[p+2] = 0;   frame[p+3] = 255; break;
			case empty:  frame[p] = 0;   frame[p+1] = 0;   frame[p+2] = 0;   frame[p+3] = 255; break;
		}
	}
}

void choose_certain_neighbor(int &i, int &j, CellType type) {
	int curr;
	switch(random_engine()%4) {
		case 0: 
			curr = (i-1+N)%N;
			if (get_celltype(matrix[curr][j]) == type) {
				i = curr;
				break;
			}
		case 1: 
			curr = (i+1)%N;
			if (get_celltype(matrix[curr][j]) == type) {
				i = curr;
				break;
			}
		case 2: 
			curr = (j-1+M)%M;
			if (get_celltype(matrix[i][curr]) == type) {
				j = curr;
				break;
			}
		case 3: 
			curr = (j+1)%M;
			if (get_celltype(matrix[i][curr]) == type) {
				j = curr;
				break;
			}
	}
}

void choose_empty_neighbor(int &neighbori, int &neighborj) { choose_certain_neighbor(neighbori, neighborj, empty); }
void choose_hering_neighbor(int &neighbori, int &neighborj) { choose_certain_neighbor(neighbori, neighborj, hering); }

chrono::duration<double> avg_thread_time;
int avg_steps_done;
int avg_herings_done;
int avg_sharks_done;
chrono::duration<double> avg_herings_time;
chrono::duration<double> avg_sharks_time;
int thread_jobs;

void iteration(size_t thread_count) {
	auto global_start = chrono::system_clock::now();
	/* cout << "*******************************" << endl; */

	//granularity of 1 (or 2, depending on how you look at it)
	int vertical_sections = thread_count;
	assert(thread_count%vertical_sections == 0);
	int horizontal_sections = thread_count / vertical_sections;
	if (vertical_sections % 2) vertical_sections *= 2;
	else horizontal_sections *= 2;
	/* if (vertical_sections < horizontal_sections) vertical_sections *= 2; */
	/* else horizontal_sections *= 2; */
	assert(N%vertical_sections == 0);
	assert(M%horizontal_sections == 0);
	int vertical_size   = N / vertical_sections;
	int horizontal_size = M / horizontal_sections;

	/* auto corners = new mutex[vertical_sections*horizontal_sections][4]; */

	auto start = chrono::system_clock::now();

	mutex summer_mutex;
	chrono::duration<double> sum_times_taken = chrono::duration<double>();

	auto iterate = [&](bool odd, int thread_idx) {
		int grid_idx = thread_idx*2 + odd;
		int row = grid_idx / horizontal_sections;
		int col = grid_idx % horizontal_sections;
		int starti = row*vertical_size;
		int startj = col*horizontal_size;
		int steps = 0;
		int hering_steps = 0;
		int shark_steps = 0;
		for (int i=starti; i<starti+vertical_size; ++i) for (int j=startj; j<startj+horizontal_size; ++j) {
			auto c = matrix[i][j];
			if (updated[i][j] || get_celltype(c) == empty) continue;
			++steps;
			bool on_up_horizontal_boundary = (i+1)%vertical_size == 0;
			bool on_down_horizontal_boundary = i%vertical_size == 0;
			bool on_left_vertical_boundary = (j+1)%horizontal_size == 0;
			bool on_right_vertical_boundary = j%horizontal_size == 0;
			bool on_corner = (on_up_horizontal_boundary || on_down_horizontal_boundary) && (on_left_vertical_boundary || on_right_vertical_boundary);

			if (on_corner) { //corner cell
				/* auto &corner = corners */
				/* 	[((i/vertical_size) % vertical_sections) * horizontal_sections + */
				/* 	(j/horizontal_size) % horizontal_sections] */
				/* 	[on_down_horizontal_boundary*2 + on_right_vertical_boundary]; */
				/* if (corner.try_lock()) corner.unlock(); */
				/* else cout << "corner already locked" << endl; */
				/* corners */
				/* 	[((i/vertical_size) % vertical_sections) * horizontal_sections + */
				/* 	(j/horizontal_size) % horizontal_sections] */
				/* 	[on_down_horizontal_boundary*2 + on_right_vertical_boundary] */
				/* 	.lock(); */
			}
			auto start = chrono::system_clock::now();
			if (get_celltype(c) == hering) {
				int neighbori = i, neighborj = j;
				choose_empty_neighbor(neighbori, neighborj);
				c = increment_counter(c, REPRODUCTION_COUNTER_MASK, HERING_REPRODUCTION_TIME);
				if (get_reproduction_counter(c) == 0) {
					matrix[i][j] = hering;
				}
				matrix[neighbori][neighborj] = c;
				updated[neighbori][neighborj] = true;
				++hering_steps;
				auto stop = chrono::system_clock::now();
				avg_herings_time += stop-start;
			} else if (get_celltype(c) == shark) {
				int neighbori = i, neighborj = j;
				choose_hering_neighbor(neighbori, neighborj);
				if (neighbori == i && neighborj == j)
					choose_empty_neighbor(neighbori, neighborj);
				c = increment_counter(c, ENERGY_COUNTER_MASK, SHARK_LIFETIME);
				c = increment_counter(c, REPRODUCTION_COUNTER_MASK, SHARK_REPRODUCTION_TIME);
				if (get_energy_counter(c) == 0) {
					matrix[i][j] = empty;
				} else {
					if (get_reproduction_counter(c) == 0) {
						matrix[i][j] = shark;
					}
					if (get_celltype(matrix[neighbori][neighborj]) == hering) {
						c = clear_energy_counter(c);
					}
					matrix[neighbori][neighborj] = c;
					updated[neighbori][neighborj] = true;
				}
				++shark_steps;
				auto stop = chrono::system_clock::now();
				avg_sharks_time += stop-start;
			}
			if (on_corner) {
				/* corners */
				/* 	[((i/vertical_size) % vertical_sections) * horizontal_sections + */
				/* 	(j/horizontal_size) % horizontal_sections] */
				/* 	[on_down_horizontal_boundary*2 + on_right_vertical_boundary] */
				/* 	.lock(); */
			}
		}
		auto stop = chrono::system_clock::now();
		summer_mutex.lock();
		sum_times_taken += stop-start;
		avg_thread_time += stop-start;
		avg_steps_done += steps;
		avg_herings_done += hering_steps;
		avg_sharks_done += shark_steps;
		++thread_jobs;
		summer_mutex.unlock();
	};


	for (int parity=0; parity<2; ++parity) {
		sum_times_taken = chrono::duration<double>();
		vector<thread> threads;
		threads.reserve(thread_count);
		start = chrono::system_clock::now();
		for (int thread_idx=0; thread_idx<thread_count; ++thread_idx) {
			threads.emplace_back(iterate, parity, thread_idx);
		}
		for (auto &th : threads) th.join();
		auto stop = chrono::system_clock::now();
		/* cout << "thread wave took time in total " << chrono::duration<double>(stop-start).count() << ", and the sum of all durations is " << sum_times_taken.count() << endl; */
	}

	/* delete[] corners; */

	fill(&updated[0][0], &updated[N][0], false);
	auto stop = chrono::system_clock::now();
	/* cout << "time for iteration " << chrono::duration<double>(stop-global_start).count() << endl; */
}

void initialize_ocean() {
	for (int i=0; i<N; ++i) {
		for (int j=0; j<M; ++j) {
			if (random_engine()%100 > 90) matrix[i][j] = random_engine()%2+1;
			else matrix[i][j] = 0;
		}
	}
}

int main(int argc, char *argv[]) {
	if (argc < 2) { cout << "You didn't specify an output gif file" << endl; return 1; }

	random_engine.seed(1);

	printf("Field is %d x %d = %d\n", N, M, N*M);

	vector<uint8_t> frame(N * M * 4, 0);

	auto fileName = argv[1];
	const int fps = 15;
	const int secs = 15;
	int delay = 100/fps;

	printf("Outputting %d frames for a duration of %d seconds which makes a total of %d frames\n", fps, secs, fps*secs);

	for (const size_t thread_count : {1,2,3,4,6,8,12}) {
		chrono::duration<double> time_worked = chrono::duration<double>();
		thread_jobs = 0;
		avg_thread_time = chrono::duration<double>();
		avg_steps_done = 0;
		avg_herings_done = 0;
		avg_sharks_done = 0;
		avg_herings_time = chrono::duration<double>();
		avg_sharks_time = chrono::duration<double>();

		/* GifWriter g; */
		/* GifBegin(&g, fileName, M, N, delay); */
		initialize_ocean();

		/* const size_t thread_count = std::thread::hardware_concurrency(); */
		/* const size_t thread_count = 1; */

		for (int i=0; i<fps*secs; ++i) {
			/* visualize_frame(frame); */
			/* GifWriteFrame(&g, frame.data(), M, N, delay); */
			auto start = chrono::system_clock::now();
			iteration(thread_count);
			auto stop = chrono::system_clock::now();
			/* cout << (double)i/(fps*secs) * 100 << "%\r"; */
			time_worked += stop - start;
		}

		avg_thread_time /= thread_jobs;
		avg_steps_done /= thread_jobs;
		avg_herings_time /= avg_herings_done;
		avg_herings_done /= thread_jobs;
		avg_sharks_time /= avg_sharks_done;
		avg_sharks_done /= thread_jobs;
		cout << "average normalized thread job time with "              << thread_count << " threads: " << avg_thread_time.count()  *thread_count << endl;
		cout << "average normalized thread job heavy steps done with "  << thread_count << " threads: " << avg_steps_done           *thread_count << endl;
		cout << "average normalized thread job hering steps done with " << thread_count << " threads: " << avg_herings_done         *thread_count << endl;
		cout << "average normalized thread job shark  steps done with " << thread_count << " threads: " << avg_sharks_done          *thread_count << endl;
		cout << "average normalized thread job hering time done with "  << thread_count << " threads: " << avg_herings_time.count() *1e8 << endl;
		cout << "average normalized thread job shark  time done with "  << thread_count << " threads: " << avg_sharks_time .count() *1e8 << endl;

		/* GifEnd(&g); */

		cout << "total time worked " << time_worked.count() << endl;
	}
}
