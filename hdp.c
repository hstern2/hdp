#include <ctype.h>
#include "hdp.h"
#include "util.h"

static int hamming_distance(const char *s1, const char *s2)
{
  int d = 0;
  while (*s1) {
    if (*s1 != *s2)
      d++;
    s1++;
    s2++;
  }
  return d;
}

static void show_match(const char *s1, const char *s2)
{
  const char *t1 = s1, *t2 = s2;
  while (*t1) {
    if (*t1 == *t2)
      printf(".");
    else
      printf("%c",*t1);
    t1++;
    t2++;
  }
  printf("\n");
  t1 = s1;
  t2 = s2;
  while (*t1) {
    if (*t1 == *t2)
      printf(".");
    else
      printf("%c", *t2);
    t1++;
    t2++;
  }
  printf("\n\n");
}

void show_matches(struct fasta_list_t f)
{
  int i, j;
  for (i = 0; i < f.n; i++)
    for (j = 0; j < f.n; j++) {
      printf("%s %s:\n", f.p[i].name, f.p[j].name);
      show_match(f.p[i].seq, f.p[j].seq);
    }
}

void show_hamming_distance_matrix(struct fasta_list_t f)
{
  int i, j;
  printf("\"\"");
  for (j = 0; j < f.n; j++)
    printf(",\"%s\"", f.p[j].name);
  printf("\n");
  for (i = 0; i < f.n; i++) {
    printf("\"\%s\"", f.p[i].name);
    for (j = 0; j < f.n; j++)
      printf(",%d", hamming_distance(f.p[i].seq, f.p[j].seq));
    printf("\n");
  }
}

#define BUFSIZE 1048576
#define MAXSUB 1024
struct fasta_list_t fasta_list_create(const char *fname, char *subset_desc)
{
  /* parse subset description */
  int nsubset = 0;
  int subset[MAXSUB][2];
  char *s;
  int seqlen = 0;
  for (s = strtok(subset_desc, ","); s && nsubset < MAXSUB; s = strtok(0, ","), nsubset++) {
    unsigned istart, iend;
    if (sscanf(s, "%u-%u", &istart, &iend) != 2)
      die("fasta_list_create: expecting <num>-<num> in subset_desc, found %s", s);
    subset[nsubset][0] = istart - 1; /* start index */
    subset[nsubset][1] = iend - istart + 1; /* length */
    if (subset[nsubset][1] < 1)
      die("fasta_list_create: end of subset must be greater than or equal to beginning of subset");
    seqlen += subset[nsubset][1];
  }
  if (nsubset >= MAXSUB)
    die("fasta_list_create: too many subsets; increase MAXSUB in hdp.c to at least %d and recompile", nsubset+1);
  /* read fasta file */
  FILE *f = safe_fopen(fname, "r");
  struct fasta_entry_t fasta[BUFSIZE];
  char buf[BUFSIZE];
  int n, nf = 0;
  while (!end_of_file(f)) {
    /* check for '>' */
    if ((buf[0] = fgetc(f)) != '>')
      die("fasta_create: %s: expecting '>', found %c", fname, buf[0]);
    /* read name */
    if (!fgets(buf, BUFSIZE, f))
      die("fasta_create: %s: unexpected end of file", fname);
    n = strlen(buf);
    if (buf[n-1] != '\n')
      die("fasta_create: %s: line too long", fname);
    buf[n-1] = '\0'; /* remove newline */
    n--;
    fasta[nf].name = safe_malloc((n+1)*sizeof(char));
    strcpy(fasta[nf].name, buf);
    /* read sequence */
    fgets(buf, BUFSIZE, f);
    n = strlen(buf);
    if (n >= BUFSIZE-1)
      die("fasta_create: %s: sequence too long; increase BUFSIZE in hdp.c to at least %d and recompile", fname, n+2);
    if (buf[n-1] == '\n') {
      buf[n-1] = '\0'; /* remove newline if present */
      n--;
    }
    fasta[nf].seq = safe_malloc((seqlen+1)*sizeof(char));
    /* copy subsets */
    int isubset, nseq = 0;
    for (isubset = 0; isubset < nsubset; isubset++) {
      const int istart = subset[isubset][0];
      const int len = subset[isubset][1];
      if (istart+len-1 >= n)
	      die("fasta_create: subset is out of range (istart: %d  len: %d  seq length: %d)", istart, len, n);
      strncpy(&fasta[nf].seq[nseq], &buf[istart], len);
      nseq += len;
    }
    fasta[nf].seq[seqlen-1] = '\0';
    nf++;
    if (nf >= BUFSIZE)
      die("fasta_create: file '%s' has too many entries; increase BUFSIZE in hdb.c to at least %d and recompile", fname, nf+1);
  }
  fclose(f);
  struct fasta_list_t fasta_list;
  fasta_list.n = nf;
  fasta_list.p = safe_malloc(nf * sizeof(struct fasta_entry_t));
  for (n = 0; n < fasta_list.n; n++)
    fasta_list.p[n] = fasta[n];
  return fasta_list;
}
#undef BUFSIZE

int fasta_list_are_lengths_identical(struct fasta_list_t fasta_list)
{
  if (fasta_list.n == 0)
    return 1;
  const int n = strlen(fasta_list.p[0].seq);
  int i;
  for (i = 1; i < fasta_list.n; i++)
    if (strlen(fasta_list.p[i].seq) != n)
      return 0;
  return 1;
}

void fasta_list_delete(struct fasta_list_t fasta_list)
{
  int i;
  for (i = 0; i < fasta_list.n; i++) {
    free(fasta_list.p[i].name);
    free(fasta_list.p[i].seq);
  }
  free(fasta_list.p);
}

void fasta_list_show(struct fasta_list_t fasta_list)
{
  int i;
  for (i = 0; i < fasta_list.n; i++)
    printf(">%s\n%s\n", fasta_list.p[i].name, fasta_list.p[i].seq);
}
