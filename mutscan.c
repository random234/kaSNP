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

struct MutScan {
  GtStrArray *vcf_arr;
  MutGene *mut_gene;  
  //~ GtFeatureNode *node, *child;
  unsigned long splice_site_interval;
};

#define GT_PHASE_CHARS \
        "012."

MutScan* mutscan_new(void) {
  MutScan *m = gt_malloc(sizeof *m);
  
  m->splice_site_interval = 0;
  
  return m;
}

GtStrArray* mutscan_get_vcf_array(MutScan *m) {
  gt_assert(m);
  return m->vcf_arr;
}

void mutscan_set_vcf_array(MutScan *m, GtStrArray *vcf_arr) {
  gt_assert(m);
  m->vcf_arr = vcf_arr;
}

MutGene* mutscan_get_mut_gene(MutScan *m) {
  gt_assert(m);
  return m->mut_gene;
}

void mutscan_set_mut_gene(MutScan *m, MutGene *g) {
  gt_assert(m);
  m->mut_gene = g;
}

unsigned long mutscan_get_splice_site_interval(MutScan *m) {
  gt_assert(m);
  return m->splice_site_interval;
}

void mutscan_set_splice_site_interval(MutScan *m, unsigned long s) {
  gt_assert(m);
  m->splice_site_interval = s;
}

unsigned long mutscan_init(MutScan *mut, GtStrArray *vcf, GtFeatureNode *fn) {
  mut->vcf_arr = vcf;
  MutGene *gene = mutgene_new();
  mut->mut_gene = gene;
   
  GtFeatureNode *node, *child;    
  GtFeatureNodeIterator *fni;
  GtFeatureNodeIterator *fni_child;
  
  GtRange rng_fn = gt_genome_node_get_range((GtGenomeNode*) fn);
  mutgene_add_content(mut->mut_gene, gt_str_new_cstr(gt_feature_node_get_type(fn)),rng_fn.start,rng_fn.end, gt_feature_node_get_phase(fn));
    
  /* get all mRNA entries in current gene */
  fni =  gt_feature_node_iterator_new_direct(fn);
  while ((node = gt_feature_node_iterator_next(fni))) {
    GtRange rng_node = gt_genome_node_get_range((GtGenomeNode*) node);
    
    /* create new child elem from current node */
    MutGene *mut_child_elem;
    mut_child_elem = mutgene_new();
    mutgene_add_content(mut_child_elem, gt_str_new_cstr(gt_feature_node_get_type(node)),rng_node.start,rng_node.end, gt_feature_node_get_phase(node));
    
    /* add new child elem to parent object */    
    mutgene_add_child(mut->mut_gene,mut_child_elem);
        
    printf("MRNACHILD: type: %s, %lu-%lu, seqid %s\n",
      gt_feature_node_get_type(node),
      rng_node.start,
      rng_node.end,
      gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)));        
    
      //~ /* get all exon & CDS entries in current gene */
      fni_child = gt_feature_node_iterator_new_direct(node);      
      while ((child = gt_feature_node_iterator_next(fni_child))) {
        GtRange rng_child = gt_genome_node_get_range((GtGenomeNode*) child);
        
        /* create new child elem for children of current node */        
        MutGene *mut_child_of_child_elem;
        mut_child_of_child_elem = mutgene_new();
        mutgene_add_content(mut_child_of_child_elem, gt_str_new_cstr(gt_feature_node_get_type(child)),rng_child.start,rng_child.end, gt_feature_node_get_phase(child));

        /* add new child elem to parent object */    
        mutgene_add_child(mut_child_elem,mut_child_of_child_elem);
        
        printf("CHILD: type: %s, %lu-%lu, seqid %s\t, phase %d\n ",
          gt_feature_node_get_type(child),
          rng_child.start,
          rng_child.end,
          gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) child)),
          gt_feature_node_get_phase(child));    
        
        //~ mutgene_delete(mut_child_of_child_elem);
      }
      //~ mutgene_delete(mut_child_elem);
  }
  
  return 0;
}


unsigned long mutscan_start_scan(MutScan *m) {
  /* check for mutations in introns */
  mutscan_intron(m);
  printf("LALA\n");
  
  return 0;
}


GtStrArray* mutscan_frms(GT_UNUSED MutScan *m) {
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr;
}

/* This function checks for nonsense & missense mutations */
GtStrArray* mutscan_miss(GT_UNUSED MutScan *m){
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr;
}

/* This function checks for mutations near splice sites */
GtStrArray* mutscan_splice(GT_UNUSED MutScan *m){
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
    
  return res_arr;
}

unsigned long mutscan_intron(MutScan *m){
  unsigned long i,j = 0;
  unsigned long var_pos = strtol(gt_str_array_get(mutscan_get_vcf_array(m),1),NULL,0);
  unsigned long ret_val;
  
  for(i = 0;i < gt_array_size(mutgene_get_children_array(mutscan_get_mut_gene(m)));i++) {
    printf("------------------------------------------------\n");
    printf("%lu\n", i);
    printf("%s \n",gt_str_get(mutgene_get_type(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i))));
    printf("%lu \n",mutgene_get_rng_start(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)));
    printf("%lu \n",mutgene_get_rng_end(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)));
    printf("%lu \n",mutgene_get_phase(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)));
    
    for(j=0;j<gt_array_size(mutgene_get_children_array(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)));j++){
      printf("%s \n",gt_str_get(mutgene_get_type(gt_array_get(mutgene_get_children_array(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)),j))));
      printf("%lu \n",mutgene_get_rng_start(gt_array_get(mutgene_get_children_array(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)),j)));
      printf("%lu \n",mutgene_get_rng_end(gt_array_get(mutgene_get_children_array(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)),j)));
      printf("%lu \n",mutgene_get_phase(gt_array_get(mutgene_get_children_array(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)),j)));
      
      if(var_pos >= mutgene_get_rng_start(gt_array_get(mutgene_get_children_array(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)),j)) && var_pos <= mutgene_get_rng_end(gt_array_get(mutgene_get_children_array(gt_array_get(mutgene_get_children_array(mutscan_get_mut_gene(m)), i)),j))) {
        ret_val = 1;
      }
      
    }
  }  
  if(!(ret_val == 1)) {
    printf("###################### We have a intron mutation at POS: %lu ###############################\n", var_pos);
  }
  return ret_val;
}


void mutscan_reset(MutScan *mut) {
  gt_assert(mut);
  gt_str_array_reset(mut->vcf_arr);
  //~ gt_genome_node_delete((GtGenomeNode*)mut->node);
}

void mutscan_delete(MutScan *mut) {
  gt_free(mut);  
}