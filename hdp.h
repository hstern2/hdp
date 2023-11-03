#ifndef HDP_H
#define HDP_H

#ifdef __cplusplus
extern "C" {
#endif

struct fasta_entry_t
{
  char *name, *seq;
};

struct fasta_list_t
{
  int n;
  struct fasta_entry_t *p;
};

struct fasta_list_t fasta_list_create(const char *fname, char *subset_desc);

void fasta_list_delete(struct fasta_list_t);
void fasta_list_show(struct fasta_list_t);
int fasta_list_are_lengths_identical(struct fasta_list_t);

void show_hamming_distance_matrix(struct fasta_list_t f);
void show_matches(struct fasta_list_t f);

#ifdef __cplusplus
}
#endif

#endif /* HDP_H */
