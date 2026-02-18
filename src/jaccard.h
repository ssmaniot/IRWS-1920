#ifndef JACCARD_H
#define JACCARD_H

void compute_jaccard(const double *v, const int *row_ptr_t, const int *col_ind_t, int no_nodes,
                     int top_K, const char *fname, const char *method);

#endif
