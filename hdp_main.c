#include "hdp.h"
#include "util.h"

int main(int argc, char *argv[])
{
  if (argc != 3)
    die("usage: %s <fasta file> <subset description>\n"
	"\n"
	"subset description should be of the form\n"
	"\n"
	"i-j,k-l,p-q\n"
	"\n"
	"e.g.\n"
	"\n"
	"%s inputFile.fasta 1-54,272-566 > out.csv"
	"\n", argv[0], argv[0]);

  struct fasta_list_t fasta_list = fasta_list_create(argv[1], argv[2]);

  if (!fasta_list_are_lengths_identical(fasta_list))
    die("hdp: %s: lengths are not identical", argv[1]);

#if 0
  fasta_list_show(fasta_list);
  show_matches(fasta_list);
#endif
  show_hamming_distance_matrix(fasta_list);

  fasta_list_delete(fasta_list);

  return 0;
}
