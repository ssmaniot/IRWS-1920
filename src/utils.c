#include "utils.h"

#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <unistd.h>

/* Helper functions */

int write_data(char path[], void *data, size_t nmemb, size_t size) {
  FILE *pdata;

  if ((pdata = fopen(path, "wb")) == NULL) {
    fprintf(stderr, " [ERROR] Cannot create file \"%s\"\n", path);
    return EXIT_FAILURE;
  }
  fwrite(data, size, nmemb, pdata);
  fclose(pdata);
  return EXIT_SUCCESS;
}

void delete_folder(char dir[]) {
  DIR *pf = opendir(dir);
  struct dirent *next_file;
  char fpath[PATH];

  while ((next_file = readdir(pf)) != NULL) {
    sprintf(fpath, "%s/%s", dir, next_file->d_name);
    remove(fpath);
  }
  closedir(pf);
  rmdir(dir);
}

void *mmap_data(char path[], size_t nmemb, size_t size) {
  int fd;
  char mmap_p[MMAP];
  void *mp;
  sprintf(mmap_p, "./%s", path);
#ifdef DEBUG
  printf("mmapping \"%s\"\n", mmap_p);
#endif
  fd = open(mmap_p, O_RDONLY);
  mp = mmap(NULL, nmemb * size, PROT_READ, MAP_SHARED, fd, 0);
  if (mp == MAP_FAILED)
    mp = NULL;
  close(fd);
  return mp;
}

void print_vec_f(double *v, int n) {
  int i;
  printf("[ ");
  for (i = 0; i < n; ++i)
    printf("%.3f ", v[i]);
  printf("]\n");
}

void print_vec_d(int *v, int n) {
  int i;
  printf("[ ");
  for (i = 0; i < n; ++i)
    printf("%d ", v[i]);
  printf("]\n");
}

void double_merge(int *from, int *to, int lo, int mid, int hi) {
  int *FL, *TL;
  int *FR, *TR;
  int i, j;
  int NL, NR;

  NL = mid - lo;
  NR = hi - mid;

  FL = (int *)malloc(sizeof(int) * NL);
  TL = (int *)malloc(sizeof(int) * NL);
  FR = (int *)malloc(sizeof(int) * NR);
  TR = (int *)malloc(sizeof(int) * NR);

  for (i = 0; i < NL; ++i) {
    FL[i] = from[lo + i];
    TL[i] = to[lo + i];
  }

  for (j = 0; j < NR; ++j) {
    FR[j] = from[mid + j];
    TR[j] = to[mid + j];
  }

  i = 0;
  j = 0;
  while (i < NL && j < NR) {
    if (TL[i] < TR[j] || (TL[i] == TR[j] && FL[i] <= FR[j])) {
      from[lo + i + j] = FL[i];
      to[lo + i + j] = TL[i];
      ++i;
    } else {
      from[lo + i + j] = FR[j];
      to[lo + i + j] = TR[j];
      ++j;
    }
  }

  while (i < NL) {
    from[lo + i + j] = FL[i];
    to[lo + i + j] = TL[i];
    ++i;
  }

  while (j < NR) {
    from[lo + i + j] = FR[j];
    to[lo + i + j] = TR[j];
    ++j;
  }

  free(FL);
  free(TL);
  free(FR);
  free(TR);
}

void double_merge_sort(int *from, int *to, int lo, int hi) {
  if (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    double_merge_sort(from, to, lo, mid);
    double_merge_sort(from, to, mid, hi);
    double_merge(from, to, lo, mid, hi);
  }
}

void sort_input_data(int *from, int *to, int n) {
  double_merge_sort(from, to, 0, n);
}

int cmp_ptr(const void *a, const void *b) {
  const double **L = (const double **)a;
  const double **R = (const double **)b;

  return (**R < **L) - (**L < **R);
}

int *index_sort_top_K(const double *v, size_t n, int top_K) {
  size_t i;
  const double **ptrs;
  int *idx;

  ptrs = (const double **)malloc(n * sizeof(const double **));
  for (i = 0; i < n; ++i)
    ptrs[i] = v + i;
  printf("Sorting pointers...");
  qsort(ptrs, n, sizeof(const double *), cmp_ptr);
  printf(" Done.\n");
  idx = (int *)malloc(sizeof(int) * top_K);
  for (i = 0; i < top_K; ++i)
    idx[i] = ptrs[n - i - 1] - v;
  free(ptrs);
  return idx;
}
