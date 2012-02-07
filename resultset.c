#include "genometools.h"
#include "resultset.h"

struct ResultSet{
  /* VCF part */
  GtStrArray *vcf_arr;
};

ResultSet* resultset_new(void) {
  ResultSet *r = gt_malloc(sizeof(ResultSet));
  return r;
}