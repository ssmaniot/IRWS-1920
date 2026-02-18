#include <stdio.h>
#include <stdlib.h>

#include "jaccard.h"
#include "utils.h"

void compute_jaccard(const double *v, const int *row_ptr_t, const int *col_ind_t, int no_nodes,
                     int top_K, const char *fname, const char *method) {
  FILE *pf = NULL;
  double **jaccard_coefficients = NULL;
  int *sorted_idx = NULL;
  int *degs = NULL;
  char fname_topk_jac[512] = {0};
  double jaccard_coefficient;
  int size_int, size_uni;
  int i, j, k;
  int ri, rj;
  int rie, rje;

  /* Creating the K x K matrixes for the top-K Jaccard Coefficients */
  jaccard_coefficients = (double **)malloc(top_K * sizeof(double *));
  for (i = 0; i < top_K; ++i) {
    jaccard_coefficients[i] = (double *)malloc(top_K * sizeof(double));
  }

  /* Computing the top-K nodes for each distribution */
  sorted_idx = index_sort_top_K(v, no_nodes, top_K);

  printf("Top-K nodes: ");
  print_vec_d(sorted_idx, top_K);

  degs = (int *)malloc(sizeof(int) * top_K);
  for (k = 0; k < top_K; ++k) {
    i = sorted_idx[k];
    degs[k] = row_ptr_t[i + 1] - row_ptr_t[i];
  }
  printf("Degree distribution: ");
  print_vec_d(degs, top_K);

  /* Creating CSV file for storing the results */
  sprintf(fname_topk_jac, "%s_%s_k%d.csv", fname, method, top_K);

  if ((pf = fopen(fname_topk_jac, "w")) == NULL) {
    fprintf(stderr, " [ERROR] cannot open output file \"%s\"\n",
            fname_topk_jac);
    exit(EXIT_FAILURE);
  }
  fprintf(pf, "n1,n2,jac\n");

  /* Computing Jaccard with a */
  for (i = 0; i < top_K; ++i) {
    for (j = i + 1; j < top_K; ++j) {
      ri = row_ptr_t[sorted_idx[i]];
      rj = row_ptr_t[sorted_idx[j]];
      rie = row_ptr_t[sorted_idx[i] + 1];
      rje = row_ptr_t[sorted_idx[j] + 1];
      size_int = 0;
      size_uni = 0;
      while (ri < rie && rj < rje) {
        if (col_ind_t[ri] < col_ind_t[rj]) {
          ++ri;
        } else if (col_ind_t[ri] > col_ind_t[rj]) {
          ++rj;
        } else {
          ++size_int;
          ++ri;
          ++rj;
        }
        ++size_uni;
      }
      while (ri < rie) {
        ++ri;
        ++size_uni;
      }
      while (rj < rje) {
        ++rj;
        ++size_uni;
      }
      jaccard_coefficient = ((double)size_int) / ((double)size_uni);
      jaccard_coefficients[i][j] = jaccard_coefficient;
      jaccard_coefficients[j][i] = jaccard_coefficient;
      printf("J(%d,%d) = %.3f\n", sorted_idx[i], sorted_idx[j],
             jaccard_coefficient);
      fprintf(pf, "%d,%d,%.3f\n", sorted_idx[i], sorted_idx[j],
              jaccard_coefficient);
    }
  }
  fclose(pf);

  free(degs);
  for (i = 0; i < top_K; ++i) {
    free(jaccard_coefficients[i]);
  }
  free(jaccard_coefficients);
  free(sorted_idx);
}
