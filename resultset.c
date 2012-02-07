#include "genometools.h"
#include "resultset.h"

struct ResultSet{
  /* VCF part */
  GtStrArray *vcf_arr;
  unsigned long var_pos;
  
  /*GFF3 part */
  GtStr *gene_name;
  
  
  
  
  /* in_exon results */
  
  
  
  /* frame shift results */
  
  /* missense nonssense results */
  
  
  
  
};

ResultSet* resultset_new(void) {
  ResultSet *r = gt_malloc(sizeof(ResultSet));
  return r;
}