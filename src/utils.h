#ifndef SHARED_H
#define SHARED_H

#include <stddef.h>

#define TOL 1.e-10
#define MAX_ITER 200
#define MOD_ITER 10
#define FNAME 256
#define DNAME 1024
#define PATH 1024
#define MMAP 2048

/* Data for compression */
typedef struct {
  int no_nodes;
  int no_edges;
  int no_danglings;
} CSR_data;

/* Helper functions */
int write_data(const char path[], const void *data, size_t nmemb, size_t size);
void delete_folder(const char dir[]);
void *mmap_data(const char path[], size_t nmemb, size_t size);
void print_vec_f(const double *v, int n);
void print_vec_d(const int *v, int n);
void double_merge(int *from, int *to, int lo, int mid, int hi);
void double_merge_sort(int *from, int *to, int lo, int hi);
void sort_input_data(int *from, int *to, int n);
int *index_sort_top_K(const double *v, size_t n, int top_K);

#endif
