#include <stdio.h>
#include <stdlib.h>
#include "genometools.h"
#include "readinput.h"

int main(int argc, char **argv) {
    const char * gff3_file;
    const char * fas_file;
    const char * vcf_file;
    
    if (argc<4)
	{
		fprintf(stderr,"Usage: %s [GFF3 file] [FASTA file] [VCF file]\n", argv[0]);
		exit(EXIT_FAILURE);
	} else {		
        gff3_file = argv[1];
        fas_file = argv[2];
        vcf_file = argv[3];
        
        printf("%s",gff3_file);
        printf("%s",fas_file);
        printf("%s",vcf_file);
        read_input(gff3_file, vcf_file, fas_file);
        
	}
    
    
    
    return 0;
}