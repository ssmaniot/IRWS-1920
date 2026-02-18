#include <dirent.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "jaccard.h"
#include "utils.h"

typedef CSR_data LCSR_data;

/* Data to save/load LCSR matrix */
char fname[FNAME] = {0};
char dir[DNAME] = {0};
char lcsr_data_p[PATH] = {0};
LCSR_data lcsr_data = {0};
int no_nodes = 0, no_edges = 0;

/*   Matrix L               Matrix L^T           */
char row_ptr_p[PATH] = {0}, row_ptr_tp[PATH] = {0};
char col_ind_p[PATH] = {0}, col_ind_tp[PATH] = {0};

/* HITS computation data */
double *a = NULL;
double *h = NULL;
char fauth[FNAME] = {0};
char fhub[FNAME] = {0};

/* LCSR matrix representation */
int *col_ind = NULL, *col_ind_t = NULL;
int *row_ptr = NULL, *row_ptr_t = NULL;

void perform_compression(const char dataset_path[FNAME]);
void compute_hits(void);

int main(int argc, char *argv[]) {
  /* File IO */
  FILE *pf;
  ssize_t bytes;

  /* Extra data */
  int err;
  int i;
  int top_K = 0;
  struct stat st = {0};

  if (argc != 2 && argc != 3) {
    fprintf(stderr,
            " [ERROR] *1* argument required: ./pagerank <arg_name> [<K>]\n");
    exit(EXIT_FAILURE);
  }

  /* Init data folder name */
  strncpy(fname, argv[1] + 5, strlen(argv[1]) - 9);
  fname[strlen(argv[1]) - 8] = '\0';
  strcpy(dir, "HITS_");
  dir[5] = '\0';
  strcat(dir, fname);
  dir[5 + strlen(fname)] = '/';
  dir[5 + strlen(fname) + 1] = '\0';

  /* Create LCSR file names */
  strcpy(row_ptr_p, dir);
  strcat(row_ptr_p, "row_ptr.bin");
  strcpy(col_ind_p, dir);
  strcat(col_ind_p, "col_ind.bin");

  /* Create transposed LCSR file names */
  strcpy(row_ptr_tp, dir);
  strcat(row_ptr_tp, "row_ptr_t.bin");
  strcpy(col_ind_tp, dir);
  strcat(col_ind_tp, "col_ind_t.bin");

  /* Create LCSR metadata file */
  strcpy(lcsr_data_p, dir);
  strcat(lcsr_data_p, "lcsr_data.bin");

  /* Create file to save HITS result */
  strcpy(fauth, fname);
  strcat(fauth, "_a.hits");
  strcpy(fhub, fname);
  strcat(fhub, "_h.hits");

  /* Check if input data has already been compressed.
   * If data has NOT yet been compressed, then perform compression */
  if (stat(dir, &st) == -1) {
    perform_compression(argv[1]);
  }

  /* Reading LCSR matrix metadata info from file */
  printf("Reading CLSR matrix data...\n");
  pf = fopen(lcsr_data_p, "rb");
  bytes = fread(&no_nodes, sizeof(lcsr_data.no_nodes), 1, pf);
  bytes = fread(&no_edges, sizeof(lcsr_data.no_edges), 1, pf);
  (void)bytes;
  fclose(pf);
  printf("no_nodes: %d\nno_edges: %d\n\n", no_nodes, no_edges);

  /* mmapping the CSR matrix data from files */
  err = 0;
  if ((row_ptr = (int *)mmap_data(row_ptr_p, sizeof(int), no_nodes + 1)) ==
      MAP_FAILED) {
    ++err;
  } else if ((row_ptr_t = (int *)mmap_data(row_ptr_tp, sizeof(int),
                                           no_nodes + 1)) == MAP_FAILED) {
    ++err;
  } else if ((col_ind = (int *)mmap_data(col_ind_p, sizeof(int), no_edges)) ==
             MAP_FAILED) {
    ++err;
  } else if ((col_ind_t = (int *)mmap_data(col_ind_tp, sizeof(int),
                                           no_edges)) == MAP_FAILED) {
    ++err;
  }

  if (err > 0) {
    fprintf(stderr, " [ERROR] data could not be mmapped from memory.\n");
    fprintf(stderr,
            "         Data is corrupted, the folder will be destroyed.\n");
    delete_folder(dir);
    /* Un-mmaping mmapped files */
    if (row_ptr == MAP_FAILED) munmap(row_ptr, (no_nodes + 1) * sizeof(int));
    if (row_ptr_t == MAP_FAILED)
      munmap(row_ptr_t, (no_nodes + 1) * sizeof(int));
    if (col_ind == MAP_FAILED) munmap(col_ind, no_edges * sizeof(int));
    if (col_ind_t == MAP_FAILED) munmap(col_ind_t, no_edges * sizeof(int));
    exit(EXIT_FAILURE);
  }

  printf("Done.\n\n");

#ifdef DEBUG
  printf("LCSR matrix\n");
  printf("---------------------\n");

  printf("col_ind: [ ");
  for (i = 0; i < no_edges; ++i) {
    printf("%d ", col_ind[i]);
  }
  printf("]\n");

  printf("row_ptr: [ ");
  for (i = 0; i < no_nodes + 1; ++i) {
    printf("%d ", row_ptr[i]);
  }
  printf("]\n\n");

  printf("Transposed LCSR matrix\n");
  printf("---------------------\n");

  printf("col_ind_t: [ ");
  for (i = 0; i < no_edges; ++i) {
    printf("%d ", col_ind_t[i]);
  }
  printf("]\n");

  printf("row_ptr_t: [ ");
  for (i = 0; i < no_nodes + 1; ++i) {
    printf("%d ", row_ptr_t[i]);
  }
  printf("]\n\n");
#endif

  /* Setting up data for HITS computation */
  a = (double *)malloc(sizeof(double) * no_nodes);
  h = (double *)malloc(sizeof(double) * no_nodes);
  for (i = 0; i < no_nodes; ++i) {
    a[i] = 1.;
    h[i] = 1.;
  }

  compute_hits();

  /* Computing top-K Jaccard coefficients */
  if (argc > 2) {
    sscanf(argv[2], "%d", &top_K);
    printf("Computing Jaccard on a\n");
    compute_jaccard(a, row_ptr_t, col_ind_t, no_nodes, top_K, fname, "a");
    printf("\nComputing Jaccard on h\n");
    compute_jaccard(h, row_ptr_t, col_ind_t, no_nodes, top_K, fname, "h");
  }

  /* un-mmapping data */
  munmap(row_ptr, (no_nodes + 1) * sizeof(int));
  munmap(row_ptr_t, (no_nodes + 1) * sizeof(int));
  munmap(col_ind, no_edges * sizeof(int));
  munmap(col_ind_t, no_edges * sizeof(int));

  /* Writing data back to memory */
  err = (write_data(fauth, (void *)a, sizeof(double), no_nodes) ==
         EXIT_FAILURE) ||
        (write_data(fhub, (void *)h, sizeof(double), no_nodes) == EXIT_FAILURE);

  /* Vectors of probability */
  free(a);
  free(h);

  /* Manage error from writing data to memory */
  if (err) {
    fprintf(stderr, " [ERROR] HITS result could not be written in memory.\n");
    if (stat(fauth, &st) == 0) remove(fauth);
    if (stat(fhub, &st) == 0) remove(fhub);
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}

void perform_compression(const char dataset_path[FNAME]) {
  /* Reading data from input file */
  FILE *pf;
  char *s = NULL;
  size_t slen = 0;
  ssize_t bytes;
  int ri, ci;
  int *from, *to;
  int f, t;
  int i;
  int err;

  printf(
      "Input file data \"%s\" is not compressed, ready to perform "
      "compression...\n\n",
      dataset_path);
  mkdir(dir, 0700);

  if ((pf = fopen(dataset_path, "r")) == NULL) {
    fprintf(stderr, " [ERROR] cannot open input file \"%s\"\n", dataset_path);
    exit(EXIT_FAILURE);
  }

  /* Parsing input file header */
  printf("Parsing input data...\n");
  bytes = getline(&s, &slen, pf);
  bytes = getline(&s, &slen, pf);
  bytes = getline(&s, &slen, pf);
  sscanf(s, "# Nodes: %d Edges: %d", &no_nodes, &no_edges);
  printf("This graph has %d nodes and %d edges\n", no_nodes, no_edges);
  bytes = getline(&s, &slen, pf);

  lcsr_data.no_nodes = no_nodes;
  lcsr_data.no_edges = no_edges;

  /* Reading data from input file */
  i = 0;
  from = (int *)malloc(sizeof(int) * no_edges);
  to = (int *)malloc(sizeof(int) * no_edges);
  while ((bytes = getline(&s, &slen, pf)) != -1) {
    sscanf(s, "%d %d", from + i, to + i);
    if (i % 10 == 0) printf("\rEdge %d/%d", i, no_edges);
    ++i;
  }
  printf("\rEdge %d/%d\n", i, no_edges);
  printf("Done\n\n");
  fclose(pf);
  free(s);

  /* LCSR matrix initialization */
  col_ind = (int *)malloc(sizeof(int) * no_edges);
  row_ptr = (int *)malloc(sizeof(int) * (no_nodes + 1));
  ri = 0;
  row_ptr[ri] = 0;
  ci = 0;

  /* Writing data in LCSR matrix */
  for (ci = 0; ci < no_edges; ++ci) {
    f = from[ci];
    t = to[ci];
    if (f > ri) {
      for (i = ri + 1; i <= f; ++i) {
        row_ptr[i] = ci;
      }
      ri = f;
    }
    col_ind[ci] = t;
  }
  while (ri < no_nodes) row_ptr[++ri] = ci;

#ifdef DEBUG
  printf("LCSR matrix\n");
  printf("---------------------\n");

  printf("col_ind: [ ");
  for (i = 0; i < no_edges; ++i) {
    printf("%d ", col_ind[i]);
  }
  printf("]\n");

  printf("row_ptr: [ ");
  for (i = 0; i < no_nodes + 1; ++i) {
    printf("%d ", row_ptr[i]);
  }
  printf("]\n\n");
#endif

  printf("Sorting edges for transposed matrix...\n");
  sort_input_data(from, to, no_edges);
  printf("Done.\n\n");

  /* Transposed LCSR matrix initialization */
  col_ind_t = (int *)malloc(sizeof(int) * no_edges);
  row_ptr_t = (int *)malloc(sizeof(int) * (no_nodes + 1));
  ri = 0;
  row_ptr_t[ri] = 0;
  ci = 0;

  /* Writing data in Transposed LCSR matrix */
  for (ci = 0; ci < no_edges; ++ci) {
    f = from[ci];
    t = to[ci];
    if (t > ri) {
      for (i = ri + 1; i <= t; ++i) {
        row_ptr_t[i] = ci;
      }
      ri = t;
    }
    col_ind_t[ci] = f;
  }
  while (ri < no_nodes) row_ptr_t[++ri] = ci;

#ifdef DEBUG
  printf("Transposed LCSR matrix\n");
  printf("---------------------\n");

  printf("col_ind_t: [ ");
  for (i = 0; i < no_edges; ++i) {
    printf("%d ", col_ind_t[i]);
  }
  printf("]\n");

  printf("row_ptr_t: [ ");
  for (i = 0; i < no_nodes + 1; ++i) {
    printf("%d ", row_ptr_t[i]);
  }
  printf("]\n\n");
#endif

  /* Writing data back to memory */
  err = (write_data(row_ptr_p, (void *)row_ptr, sizeof(int), no_nodes + 1) ==
         EXIT_FAILURE) ||
        (write_data(col_ind_p, (void *)col_ind, sizeof(int), no_edges) ==
         EXIT_FAILURE) ||
        (write_data(row_ptr_tp, (void *)row_ptr_t, sizeof(int), no_nodes + 1) ==
         EXIT_FAILURE) ||
        (write_data(col_ind_tp, (void *)col_ind_t, sizeof(int), no_edges) ==
         EXIT_FAILURE) ||
        (write_data(lcsr_data_p, (void *)&lcsr_data, sizeof(LCSR_data), 1) ==
         EXIT_FAILURE);

  /* Input data */
  free(from);
  free(to);
  from = NULL;
  to = NULL;
  /* CSR data structure */
  free(col_ind);
  free(row_ptr);
  free(col_ind_t);
  free(row_ptr_t);
  col_ind = NULL;
  row_ptr = NULL;
  col_ind_t = NULL;
  row_ptr_t = NULL;

  /* Manage error from writing data to memory */
  if (err) {
    delete_folder(dir);
    fprintf(stderr, " [ERROR] data could not be written in memory.\n");
    exit(EXIT_FAILURE);
  }
}

void compute_hits(void) {
  /* HITS computation */
  double a_dist = DBL_MAX, h_dist = DBL_MAX;
  double sum;
  int iter = 0;
  int ri, ci;
  int i;
  double *a_new = (double *)malloc(sizeof(double) * no_nodes);
  double *h_new = (double *)malloc(sizeof(double) * no_nodes);

  /* Time elapsed data */
  clock_t begin, end;
  double elapsed_time;

  /* Computing HITS */
  printf("Computing HITS...\n");
  begin = clock();
  while ((a_dist > TOL || h_dist > TOL) && iter < MAX_ITER) {
    if (iter % MOD_ITER == 0) {
      printf("\riter %d", iter);
#ifdef DEBUG
      printf("\n");
      printf("a: ");
      print_vec_f(a, no_nodes);
      printf("h: ");
      print_vec_f(h, no_nodes);
#endif
    }

    /* a_new = Lt @ h, h_new = L @ a */
    for (ri = 0; ri < no_nodes; ++ri) {
      a_new[ri] = 0.;
      for (ci = row_ptr_t[ri]; ci < row_ptr_t[ri + 1]; ++ci) {
        a_new[ri] += h[col_ind_t[ci]];
      }
      h_new[ri] = .0;
      for (ci = row_ptr[ri]; ci < row_ptr[ri + 1]; ++ci) {
        h_new[ri] += a[col_ind[ci]];
      }
    }

    /* Normalization step */
    sum = 0.;
    for (i = 0; i < no_nodes; ++i) sum += a_new[i];
    for (i = 0; i < no_nodes; ++i) a_new[i] /= sum;
    sum = 0.;
    for (i = 0; i < no_nodes; ++i) sum += h_new[i];
    for (i = 0; i < no_nodes; ++i) h_new[i] /= sum;

    /* Computing distance between current and old a/h */
    a_dist = 0.;
    h_dist = 0.;
    for (i = 0; i < no_nodes; ++i) {
      a_dist += (a[i] - a_new[i]) * (a[i] - a_new[i]);
      h_dist += (h[i] - h_new[i]) * (h[i] - h_new[i]);
    }
    a_dist = sqrt(a_dist);
    h_dist = sqrt(h_dist);

    /* Copy new values in a/h */
    for (i = 0; i < no_nodes; ++i) {
      a[i] = a_new[i];
      h[i] = h_new[i];
    }

    ++iter;
  }
  end = clock();
  printf("\riter %d\n", iter);
#ifdef DEBUG
  printf("a: ");
  print_vec_f(a, no_nodes);
  printf("h: ");
  print_vec_f(h, no_nodes);
#endif
  printf("Done.\n\n");

  printf("Proof of correctness:\n");
  sum = 0.;
  for (i = 0; i < no_nodes; ++i) {
    sum += a[i];
  }
  printf("sum(a) = %f\n", sum);
  sum = 0.;
  for (i = 0; i < no_nodes; ++i) {
    sum += h[i];
  }
  printf("sum(h) = %f\n\n", sum);

  elapsed_time = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Elapsed time: %.3fs\n", elapsed_time);

  free(a_new);
  free(h_new);
}
