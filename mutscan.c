#include "genometools.h"
#include "core/str_array_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "extended/genome_node_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "mutscan.h"

struct MutScan {
  GtStrArray *vcf_arr;
  GtFeatureNode *node, *child;
  
};

MutScan* mutscan_new(void) {
  MutScan *mut = gt_malloc(sizeof (MutScan));
  
  return mut;
}

GtStrArray* mutscan_frms(MutScan *mut, GtStrArray *vcf, GtFeatureNode *fn) {
  mut->vcf_arr = vcf;
  mut->node = fn;
  
  GtFeatureNodeIterator *fni;
  fni = gt_feature_node_iterator_new(mut->node);  
  
  fni =  gt_feature_node_iterator_new_direct(mut->node);      
  while ((mut->child = gt_feature_node_iterator_next(fni))) {
    GtRange rng_child = gt_genome_node_get_range((GtGenomeNode*) mut->child);
    printf("CHILD: type: %s, %lu-%lu, seqid %s\n",
      gt_feature_node_get_type(mut->child),
      rng_child.start,
      rng_child.end,
      gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) mut->child)));        
  }       
  return vcf;
}

GtStrArray* mutscan_miss(MutScan *mut){
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr
}

GtStrArray* mutscan_splice(MutScan *mut){
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr
}




void mutscan_reset(MutScan *mut) {
  gt_assert(mut);
  gt_str_array_reset(mut->vcf_arr);
  gt_genome_node_delete((GtGenomeNode*)mut->node);
}

void mutscan_delete(MutScan *mut) {
  gt_free(mut);  
}