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
#include "core/splitter_api.h"
#include "core/phase_api.h"
#include "core/translator_api.h"
#include "core/codon_iterator_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/trans_table_api.h"
#include "mutscan.h"
#include "mutgene.h"
#include "resultset.h"

struct MutScan {
  GtStrArray *vcf_arr;
  unsigned long file_chromosome; /* chromosome file mapping for desc file */
  MutGene *mut_gene;  
  GtEncseq *encseq;
  GtEncseqReader *encseq_read;
  GtError *err;
  unsigned long splice_site_interval;
};

#define GT_PHASE_CHARS \
        "012."

MutScan* mutscan_new(void) {
  MutScan *m = gt_malloc(sizeof *m);
  
  m->splice_site_interval = 0;
  m->err = gt_error_new();
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

unsigned long mutscan_init(MutScan *mut, GtStrArray *vcf, unsigned long file_chrom, GtFeatureNode *fn, GtEncseq *en, GtEncseqReader *rd) {
  mut->vcf_arr = vcf;
  mut->file_chromosome = file_chrom;
  MutGene *gene = mutgene_new();
  mut->mut_gene = gene;
  mut->encseq = en;
  mut->encseq_read = rd;
  
     
  GtFeatureNode *node, *child;    
  GtFeatureNodeIterator *fni;
  GtFeatureNodeIterator *fni_child;
  
  GtRange rng_fn = gt_genome_node_get_range((GtGenomeNode*) fn);
  
  mutgene_add_content(mut->mut_gene, gt_str_new_cstr(gt_feature_node_get_attribute(fn,"id")), gt_str_new_cstr(gt_feature_node_get_attribute(fn,"name")) ,gt_str_new_cstr(gt_feature_node_get_type(fn)),rng_fn.start,rng_fn.end, gt_feature_node_get_phase(fn));
    
  /* get all mRNA entries in current gene */
  fni =  gt_feature_node_iterator_new_direct(fn);
  while ((node = gt_feature_node_iterator_next(fni))) {
    GtRange rng_node = gt_genome_node_get_range((GtGenomeNode*) node);
    
    /* create new child elem from current node */
    MutGene *mut_child_elem;
    mut_child_elem = mutgene_new();
    
    //~ printf("\t Genename:%s\n\n", gt_feature_node_get_attribute(node,"Name"));
    mutgene_add_content(mut_child_elem, gt_str_new_cstr(gt_feature_node_get_attribute(node,"ID")), gt_str_new_cstr(gt_feature_node_get_attribute(node,"Name")) , gt_str_new_cstr(gt_feature_node_get_type(node)),rng_node.start,rng_node.end, gt_feature_node_get_phase(node));
    
    /* add new child elem to parent object */    
    mutgene_add_child(mut->mut_gene,mut_child_elem);
        
    //~ printf("MRNACHILD: type: %s, %lu-%lu, seqid %s\n",
      //~ gt_feature_node_get_type(node),
      //~ rng_node.start,
      //~ rng_node.end,
      //~ gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)));        
    
      //~ /* get all exon & CDS entries in current gene */
      fni_child = gt_feature_node_iterator_new_direct(node);      
      while ((child = gt_feature_node_iterator_next(fni_child))) {
        GtRange rng_child = gt_genome_node_get_range((GtGenomeNode*) child);
        
        /* create new child elem for children of current node */        
        MutGene *mut_child_of_child_elem;
        mut_child_of_child_elem = mutgene_new();
        mutgene_add_content(mut_child_of_child_elem, gt_str_new_cstr(gt_feature_node_get_attribute(child,"ID")), gt_str_new_cstr(gt_feature_node_get_attribute(child,"Name")), gt_str_new_cstr(gt_feature_node_get_type(child)),rng_child.start,rng_child.end, gt_feature_node_get_phase(child));

        /* add new child elem to parent object */    
        mutgene_add_child(mut_child_elem,mut_child_of_child_elem);
        
        //~ printf("CHILD: type: %s, %lu-%lu, seqid %s\t, phase %d\n ",
          //~ gt_feature_node_get_type(child),
          //~ rng_child.start,
          //~ rng_child.end,
          //~ gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) child)),
          //~ gt_feature_node_get_phase(child));    
        
        //~ mutgene_delete(mut_child_of_child_elem);
      }
      //~ mutgene_delete(mut_child_elem);
  }
  
  return 0;
}


ResultSet* mutscan_start_scan(MutScan *m) {
  //~ int i = 0;
  ResultSet *r = resultset_new();
  
  /* check for mutations in introns */
  /* maybe checking for mutations in exons may be more useful as one could use the information to stop some function calls of subsequent analysis */  
  mutscan_exon(m,r);
  //~ for(i=0;i<gt_str_array_size(exon_res);i++) {
    //~ printf("%s \t",gt_str_array_get(exon_res, i));
  //~ }
  //~ printf("\n");
  
  /* check for mutations in frames */
  mutscan_frame(m, r);
  //~ for(i=0;i<gt_str_array_size(frame_res);i++) {
    //~ printf("%s \t",gt_str_array_get(frame_res, i));
  //~ }
  mutscan_miss(m, r);
  
  printf("\n");
  return 0;
}

unsigned long mutscan_exon(MutScan *m,  ResultSet *r){
  unsigned long i,j,had_err = 0;
  unsigned long var_pos = strtol(gt_str_array_get(mutscan_get_vcf_array(m),1),NULL,0);
  resultset_set_var_pos(r,var_pos);
  resultset_set_gene_name(r, mutgene_get_gene_name(mutscan_get_mut_gene(m)));
  
  printf("------- mutscan_exon() -------\n");
  //~ printf("ID: %s\n",gt_str_get(mutgene_get_id(mutscan_get_mut_gene(m))));
  
  GtArray *mrna_arr = mutgene_get_children_array(mutscan_get_mut_gene(m));
  for(i = 0;i < gt_array_size(mutgene_get_children_array(mutscan_get_mut_gene(m)));i++) {
    MutGene *mrna_elem = gt_array_get(mrna_arr, i);
        printf("ID: %s\t",gt_str_get(mutgene_get_id(mrna_elem)));
        printf("Genename: %s\t",gt_str_get(mutgene_get_gene_name(mrna_elem)));
        printf("Type: %s\n",gt_str_get(mutgene_get_type(mrna_elem)));
    
    GtArray *mrna_child_arr = mutgene_get_children_array(mrna_elem);
    for(j=0;j<gt_array_size(mrna_child_arr);j++){
      MutGene *mrna_child_elem = gt_array_get(mrna_child_arr, j);
            printf("%s \t",gt_str_get(mutgene_get_type(mrna_child_elem)));
      printf("%lu \t",mutgene_get_rng_start(mrna_child_elem));
      printf("%lu \t",mutgene_get_rng_end(mrna_child_elem));
      printf("%lu \t",mutgene_get_phase(mrna_child_elem));
        printf("ID: %s\t",gt_str_get(mutgene_get_id(mrna_child_elem)));
        printf("Genename: %s\t",gt_str_get(mutgene_get_gene_name(mrna_child_elem)));
      printf("Type: %s\n",gt_str_get(mutgene_get_type(mrna_child_elem)));
      
      if(var_pos >= mutgene_get_rng_start(mrna_child_elem) && var_pos <= mutgene_get_rng_end(mrna_child_elem)) {
        /* adding mrna_id to resultset to enable subsequent functions to skip non-exonic variations */

        
        
        resultset_add_mrna_id(r, mutgene_get_id(mrna_elem));        
      } 
      mrna_child_elem = NULL;
    }
    mrna_elem = NULL;
  }
  return had_err;
}

unsigned long mutscan_frame(MutScan *m, ResultSet *r) {
  unsigned long i,j,k,had_err = 0;
  resultset_set_vcf_array(r,mutscan_get_vcf_array(m));
  GtStrArray *ref;
  GtStr *ref_temp;
  GtSplitter *ref_split;
  GtStrArray *alt;
  GtStr *alt_temp;
  GtSplitter *alt_split;
  
  ref = gt_str_array_new();
  ref_temp = gt_str_new();
  alt_temp = gt_str_new();
  alt = gt_str_array_new();
  
  ref_split = gt_splitter_new();
  alt_split = gt_splitter_new();
  
  ref_temp = gt_str_new_cstr(gt_str_array_get(resultset_get_vcf_array(r),3));
  alt_temp = gt_str_new_cstr(gt_str_array_get(resultset_get_vcf_array(r),4));
  
  gt_splitter_split(ref_split, gt_str_get(ref_temp),
                    gt_str_length(ref_temp), ',');
  gt_splitter_split(alt_split, gt_str_get(alt_temp),
                    gt_str_length(alt_temp), ',');

  for(j=0;j<gt_splitter_size(ref_split);j++) {
    gt_str_array_add_cstr(ref,gt_splitter_get_token(ref_split,j));
  }
  for(j=0;j<gt_splitter_size(alt_split);j++) {
    gt_str_array_add_cstr(alt,gt_splitter_get_token(alt_split,j));
  }  
  
  printf("ref %s alt %s", gt_str_get(ref_temp), gt_str_get(alt_temp));
  gt_str_reset(ref_temp);
  gt_str_reset(alt_temp);
  
  printf("------- mutscan_frame() -------\n");
  GtArray *mrna_arr = mutgene_get_children_array(mutscan_get_mut_gene(m));
  for(i = 0;i < gt_array_size(mutgene_get_children_array(mutscan_get_mut_gene(m)));i++) {
    MutGene *mrna_elem = gt_array_get(mrna_arr, i);
    
      //~ printf("%s \t",gt_str_get(mutgene_get_id(mrna_elem)));
      //~ printf("%lu \t",mutgene_get_rng_start(mrna_elem));
      //~ printf("%lu \t",mutgene_get_rng_end(mrna_elem));
      //~ printf("%lu \t",mutgene_get_phase(mrna_elem));
    
    if(resultset_check_mrna_ids(r, mutgene_get_id(mrna_elem)) == 0) {
      printf("found a mRNA with Variation id: %s at pos: %lu\n", gt_str_get(mutgene_get_id(mrna_elem)), resultset_get_var_pos(r));
      
      for(j=0;j<gt_str_array_size(ref);j++) {
        gt_str_set(ref_temp,gt_str_array_get(ref,j));
        
        for(k=0;k<gt_str_array_size(alt);k++) {
          gt_str_set(alt_temp,gt_str_array_get(alt,k));
          if((gt_str_length(ref_temp) - gt_str_length(alt_temp)) == 0 || (gt_str_length(ref_temp) - gt_str_length(alt_temp)) % 3 == 0) {
            printf("found frameshifting variation in gene %s \n", gt_str_get(mutgene_get_gene_name(mrna_elem)));
          }
          gt_str_reset(alt_temp);
        }
        gt_str_reset(ref_temp);
      }
    }
    mrna_elem = NULL;
  }
  return had_err;
}

/* This function checks for nonsense & missense mutations */
unsigned long mutscan_miss(MutScan *m,  ResultSet *r){
  unsigned long i,j,k,l,range,rel_var_pos,had_err = 0;
  char amino;  
  char amino_mut;
  GtStr *nucl_mut;
  GtStr *temp_arr; 
  GtSplitter *alt_split;
  GtEncseqReader *read;
  GtStr *prot_seq;
  GtTransTable *trans_t;
  read = m->encseq_read;  
  prot_seq = gt_str_new();
  alt_split = gt_splitter_new();
  
  trans_t = gt_trans_table_new_standard(m->err);  
  
  printf("------- mutscan_miss() -------\n");  
  GtArray *mrna_arr = mutgene_get_children_array(mutscan_get_mut_gene(m));
  for(i = 0;i < gt_array_size(mutgene_get_children_array(mutscan_get_mut_gene(m)));i++) {
    MutGene *mrna_elem = gt_array_get(mrna_arr, i);    
    GtArray *mrna_child_arr = mutgene_get_children_array(mrna_elem);
    for(j=0;j<gt_array_size(mrna_child_arr);j++){
      MutGene *mrna_child_elem = gt_array_get(mrna_child_arr, j);
      gt_encseq_reader_reinit_with_readmode(read,m->encseq,0,gt_encseq_seqstartpos(m->encseq,m->file_chromosome) + mutgene_get_rng_start(mrna_child_elem)-1);
      
  
      if(resultset_get_var_pos(r) < mutgene_get_rng_end(mrna_child_elem) && resultset_get_var_pos(r) > mutgene_get_rng_start(mrna_child_elem)) {
        range = mutgene_get_rng_end(mrna_child_elem) - mutgene_get_rng_start(mrna_child_elem);
        printf("Start of exon: %lu Range: %lu\n",mutgene_get_rng_start(mrna_child_elem),range);
      
        for(k=0;k<range;k++){
          gt_str_append_char(prot_seq,gt_encseq_reader_next_decoded_char(read));
        }
        resultset_add_prot_seq(r, prot_seq);

        rel_var_pos = resultset_get_var_pos(r) - mutgene_get_rng_start(mrna_child_elem) + mutgene_get_phase(mrna_child_elem);
        printf("Variation pos: %lu\n",resultset_get_var_pos(r));
        printf("Relative Variation pos: %lu\n",rel_var_pos);
        if(gt_trans_table_translate_codon(trans_t,gt_str_get(prot_seq)[rel_var_pos], gt_str_get(prot_seq)[rel_var_pos+1], gt_str_get(prot_seq)[rel_var_pos+2], &amino, m->err)== 0) {
          printf("ORIGINAL AMINOACID: %c \n",amino);        
        } else {
          gt_error_set(m->err,"error during ORIGINAL AA translation");
        }
        
        temp_arr = gt_str_new_cstr(gt_str_array_get(resultset_get_vcf_array(r),4));        
        gt_splitter_split(alt_split, gt_str_get(temp_arr),gt_str_length(temp_arr), ',');
        printf("splitter size %lu: \n",gt_splitter_size(alt_split));
        for(l=0;l<gt_splitter_size(alt_split);l++) {
          nucl_mut = gt_str_new_cstr(gt_splitter_get_token(alt_split, l));          
          if(gt_str_length(nucl_mut)==1) {
            printf("added char %c: \n",gt_str_get(prot_seq)[rel_var_pos+1]);        
            gt_str_append_char(nucl_mut,gt_str_get(prot_seq)[rel_var_pos+1]);
          }
          if(gt_str_length(nucl_mut)==2) {
            printf("added char %c: \n",gt_str_get(prot_seq)[rel_var_pos+2]);        
            gt_str_append_char(nucl_mut,gt_str_get(prot_seq)[rel_var_pos+2]);
          }
          printf("nucl_mut: %s \n",gt_str_get(nucl_mut));        
          if(gt_trans_table_translate_codon(trans_t,gt_str_get(nucl_mut)[0], gt_str_get(nucl_mut)[1], gt_str_get(nucl_mut)[2], &amino_mut, m->err)== 0) {
            printf("MUTATED AMINOACID: %c \n",amino_mut);        
          } else {            
            gt_error_set(m->err,"error during ALTERNATE AA translation");
          }
          gt_str_reset(nucl_mut);
        }
        gt_str_reset(temp_arr);
        gt_splitter_reset(alt_split);
        gt_str_reset(prot_seq);
        mrna_child_elem = NULL;
      }
    }
    
    mrna_elem = NULL;
  }
  return had_err;
}

/* This function checks for mutations near splice sites */
unsigned long mutscan_splice(MutScan *m, GT_UNUSED ResultSet *r){
  unsigned long i,j,had_err = 0;  
  printf("------- mutscan_splice() -------\n");
  GtArray *mrna_arr = mutgene_get_children_array(mutscan_get_mut_gene(m));
  for(i = 0;i < gt_array_size(mutgene_get_children_array(mutscan_get_mut_gene(m)));i++) {
    MutGene *mrna_elem = gt_array_get(mrna_arr, i);
    
    GtArray *mrna_child_arr = mutgene_get_children_array(mrna_elem);
    GtStr *res_str = gt_str_new();
    for(j=0;j<gt_array_size(mrna_child_arr);j++){
      GT_UNUSED MutGene *mrna_child_elem = gt_array_get(mrna_child_arr, j);
    
      mrna_child_elem = NULL;
    }
    gt_str_delete(res_str);
    mrna_elem = NULL;
  }
  return had_err;
}




void mutscan_reset(MutScan *mut) {
  gt_assert(mut);
  gt_str_array_reset(mut->vcf_arr);
  //~ gt_genome_node_delete((GtGenomeNode*)mut->node);
}

void mutscan_delete(MutScan *mut) {
  gt_free(mut);  
}