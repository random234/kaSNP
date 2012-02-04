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
  //~ GtFeatureNode *node, *child;
  
  
  unsigned long splice_site_interval;
};

MutScan* mutscan_new(void) {
  MutScan *mut = gt_malloc(sizeof (MutScan));
  
  return mut;
}

unsigned long mutscan_init(MutScan *mut, GtStrArray *vcf, GtFeatureNode *fn) {
  mut->vcf_arr = vcf;
  GtFeatureNode *node, *child;    
  GtFeatureNodeIterator *fni;
  GtFeatureNodeIterator *fni_child;
  
  fni =  gt_feature_node_iterator_new_direct(fn);      
  while ((node = gt_feature_node_iterator_next(fni))) {
    GtRange rng_node = gt_genome_node_get_range((GtGenomeNode*) node);
    printf("CHILD: type: %s, %lu-%lu, seqid %s\n",
      gt_feature_node_get_type(node),
      rng_node.start,
      rng_node.end,
      gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)));        
      fni_child = gt_feature_node_iterator_new(node);
      while ((child = gt_feature_node_iterator_next(fni_child))) {
        GtRange rng_child = gt_genome_node_get_range((GtGenomeNode*) child);
        printf("CHILD: type: %s, %lu-%lu, seqid %s\n",
          gt_feature_node_get_type(child),
          rng_child.start,
          rng_child.end,
          gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) child)));    
      }
  }
  return 0;
}


GtStrArray* mutscan_frms(GT_UNUSED MutScan *mut) {
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr;
}

/* This function checks for nonsense & missense mutations */
GtStrArray* mutscan_miss(GT_UNUSED MutScan *mut){
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr;
}

/* This function checks for mutations near splice sites */
GtStrArray* mutscan_splice(GT_UNUSED MutScan *mut){
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr;
}

GtStrArray* mutscan_intron(GT_UNUSED MutScan *mut){
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
  
  
  
  
  return res_arr;
}


void mutscan_reset(MutScan *mut) {
  gt_assert(mut);
  gt_str_array_reset(mut->vcf_arr);
  //~ gt_genome_node_delete((GtGenomeNode*)mut->node);
}

void mutscan_delete(MutScan *mut) {
  gt_free(mut);  
}