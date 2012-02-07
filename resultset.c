#include "genometools.h"
#include "resultset.h"

struct ResultSet{
  /* VCF part */
  GtStrArray *vcf_arr;
  unsigned long var_pos;
  
  /*GFF3 part */
  GtStr *gene_name;
  GtStr *id;
  
  /* in_exon results */
    
  /* frame shift results */
  
  /* missense nonssense results */
  
  
  
  
};

ResultSet* resultset_new(void) {
  ResultSet *r = gt_malloc(sizeof(ResultSet));
  return r;
}

void resultset_set_vcf_array(ResultSet *r, GtStrArray *v) {
  gt_assert(r);
  r->vcf_arr = v;
}

GtStrArray * resultset_get_vcf_array(ResultSet *r) {
  gt_assert(r);
  return r->vcf_arr;
}

void resultset_set_var_pos(ResultSet *r, unsigned long p) {
  gt_assert(r);
  r->var_pos = p;
}

unsigned long resultset_get_var_pos(ResultSet *r) {
  gt_assert(r);
  return r->var_pos;
}

void resultset_set_id(ResultSet *r, GtStr *i) {
  gt_assert(r);
  r->id = i;
}

GtStr * resultset_get_id(ResultSet *r) {
  gt_assert(r);
  return r->id;
}

void resultset_set_gene_name(ResultSet *r, GtStr *gn) {
  gt_assert(r);
  r->gene_name = gn;
}

GtStr * resultset_get_gene_name(ResultSet *r) {
  gt_assert(r);
  return r->gene_name;
}
