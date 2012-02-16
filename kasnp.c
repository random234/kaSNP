#include <stdio.h>
#include <stdlib.h>
#include "genometools.h"
#include "readinput.h"

int main(int argc, char **argv) {
    const char *gff3_file;
    const char *encseq_file;
    const char *vcf_file;
    const char *out_file;
    unsigned long splice_site_range;
    
    if (argc<5)
	{
		fprintf(stderr,"Usage: %s [GFF3 file] [FASTA file] [VCF file] [OUTPUT VCF file] [Splice Site interval]\n", argv[0]);
		exit(EXIT_FAILURE);
	} else {		
        gff3_file = argv[1];
        encseq_file = argv[2];
        vcf_file = argv[3];
        out_file = argv[4];
        
        splice_site_range = strtol(argv[5], NULL,0);

        read_input(gff3_file, vcf_file, encseq_file, out_file, splice_site_range);
	}
    return 0;
}