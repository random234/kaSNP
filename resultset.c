#include <string.h>
#include "genometools.h"
#include "core/str_array_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "extended/genome_node_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/phase_api.h"
#include "mutscan.h"
#include "mutgene.h"
#include "resultset.h"

struct ResultSet{
  /* VCF part */
  GtStrArray *vcf_arr;
  unsigned long var_pos;
  
  /*GFF3 part */
  GtStr *gene_name;
  GtStr *id;
  
  /* in_exon results */
  GtArray *mrna_ids;
  //~ GtStrArray *mrna_ids;    
  
  /* frame shift results */
  unsigned long frms;
  /* missense nonssense results */

};

ResultSet* resultset_new(void) {
  ResultSet *r = gt_malloc(sizeof(*r));
  r->mrna_ids = gt_array_new(sizeof (GtStr *));
  //~ r->mrna_ids = gt_str_array_new();
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

void resultset_add_mrna_id(ResultSet *r, GtStr *mi) {
  gt_assert(r);
  gt_array_add(r->mrna_ids, mi);
  //~ gt_str_array_add(r->mrna_ids, mi);
}

GtArray * resultset_get_mrna_ids(ResultSet *r) {
  gt_assert(r);
  return r->mrna_ids;
}

unsigned long resultset_check_mrna_ids(ResultSet *r, GtStr *id) {
  int i;
  unsigned long res = 1;
  
  for(i = 0;i < gt_array_size(r->mrna_ids); i++) {
    if(strcmp(gt_array_get(r->mrna_ids, 0), gt_str_get(id)) == 0) {
      res = 0;
    }    
  }  
  return res;
}


