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
  int i = 0;
  /* check for mutations in introns */
  /* maybe checking for mutations in exons may be more useful as one could use the information to stop some function calls of subsequent analysis */  
  GtStrArray *intron_res = mutscan_exon(m);
  for(i=0;i<gt_str_array_size(intron_res);i++) {
    printf("%s \t",gt_str_array_get(intron_res, i));
  }
  printf("\n");
  
  /* check for mutations in frames */
  //~ GtStrArray *frame_res = mutscan_frame(m);
  for(i=0;i<gt_str_array_size(intron_res);i++) {
    //~ printf("%s \t",gt_str_array_get(frame_res, i));
  }
  printf("\n");
  return 0;
}


GtStrArray* mutscan_frame(GT_UNUSED MutScan *m) {
  unsigned long i,j = 0;
  unsigned long var_pos = strtol(gt_str_array_get(mutscan_get_vcf_array(m),1),NULL,0);
  
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
  
  printf("------- mutscan_frms() -------\n");
  GtArray *mrna_arr = mutgene_get_children_array(mutscan_get_mut_gene(m));
  for(i = 0;i < gt_array_size(mutgene_get_children_array(mutscan_get_mut_gene(m)));i++) {
    MutGene *mrna_elem = gt_array_get(mrna_arr, i);
    //~ printf("%lu\n", i);
    //~ printf("%s \n",gt_str_get(mutgene_get_type(mrna_elem)));
    //~ printf("%lu \n",mutgene_get_rng_start(mrna_elem));
    //~ printf("%lu \n",mutgene_get_rng_end(mrna_elem));
    //~ printf("%lu \n",mutgene_get_phase(mrna_elem));
    
    GtArray *mrna_child_arr = mutgene_get_children_array(mrna_elem);
    for(j=0;j<gt_array_size(mrna_child_arr);j++){
      MutGene *mrna_child_elem = gt_array_get(mrna_child_arr, j);
      //~ printf("%s \n",gt_str_get(mutgene_get_type(mrna_child_elem)));
      //~ printf("%lu \n",mutgene_get_rng_start(mrna_child_elem));
      //~ printf("%lu \n",mutgene_get_rng_end(mrna_child_elem));
      //~ printf("%lu \n",mutgene_get_phase(mrna_child_elem));
    
      if(var_pos >= mutgene_get_rng_start(mrna_child_elem) && var_pos <= mutgene_get_rng_end(mrna_child_elem)) {
        
      }
      mrna_child_elem = NULL;
    }
    mrna_elem = NULL;
  }
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

GtStrArray* mutscan_exon(MutScan *m){
  unsigned long i,j = 0;
  unsigned long var_pos = strtol(gt_str_array_get(mutscan_get_vcf_array(m),1),NULL,0);
  GtStrArray *res_arr;
  res_arr = gt_str_array_new();
  GtStr *temp = gt_str_new();
  gt_str_append_ulong(temp, var_pos);
  
  gt_str_array_add_cstr(res_arr,"in_exon");
  gt_str_array_add(res_arr, temp);
  
  printf("------- mutscan_exon() -------\n");
  GtArray *mrna_arr = mutgene_get_children_array(mutscan_get_mut_gene(m));
  for(i = 0;i < gt_array_size(mutgene_get_children_array(mutscan_get_mut_gene(m)));i++) {
    MutGene *mrna_elem = gt_array_get(mrna_arr, i);
    
    
    //~ printf("%lu\n", i);
    //~ printf("%s \n",gt_str_get(mutgene_get_type(mrna_elem)));
    //~ printf("%lu \n",mutgene_get_rng_start(mrna_elem));
    //~ printf("%lu \n",mutgene_get_rng_end(mrna_elem));
    //~ printf("%lu \n",mutgene_get_phase(mrna_elem));
    
    GtArray *mrna_child_arr = mutgene_get_children_array(mrna_elem);
    GtStr *res_str = gt_str_new();
    for(j=0;j<gt_array_size(mrna_child_arr);j++){
      MutGene *mrna_child_elem = gt_array_get(mrna_child_arr, j);
      //~ printf("%s \n",gt_str_get(mutgene_get_type(mrna_child_elem)));
      //~ printf("%lu \n",mutgene_get_rng_start(mrna_child_elem));
      //~ printf("%lu \n",mutgene_get_rng_end(mrna_child_elem));
      //~ printf("%lu \n",mutgene_get_phase(mrna_child_elem));
    
      if(var_pos >= mutgene_get_rng_start(mrna_child_elem) && var_pos <= mutgene_get_rng_end(mrna_child_elem)) {
        gt_str_append_str(res_str, mutgene_get_type(mrna_elem));
        gt_str_append_ulong(res_str,i+1);
        gt_str_array_add(res_arr, res_str);
      } 
      mrna_child_elem = NULL;
    }
    gt_str_delete(res_str);
    mrna_elem = NULL;
  }
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