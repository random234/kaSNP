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
#include "core/array_api.h"
#include "mutgene.h"


struct MutGene {  
  GtStr *id;
  GtStr *gene_name;
  GtStr *type;  
  unsigned long rng_start;
  unsigned long rng_end;
  unsigned long phase;
  /* may not be needed as we could also get the array size from GtArray */
  unsigned long child_size;
  GtArray *children;  
};


MutGene* mutgene_new(void) {
  MutGene *g = gt_malloc(sizeof *g);
  g->id = gt_str_new();
  g->gene_name = gt_str_new();
  g->type = gt_str_new();
  g->rng_start = 0;
  g->rng_end = 0;
  g->phase = 0;
  g->children = gt_array_new(sizeof (MutGene));
  g->child_size = 0;
  return g;
}

GtStr* mutgene_get_id(MutGene *g) {
  gt_assert(g);
  return g->id;  
}

void mutgene_set_id(MutGene *g, GtStr *i) {
  gt_assert(g);
  g->id = i;  
}

GtStr* mutgene_get_gene_name(MutGene *g) {
  gt_assert(g);
  return g->gene_name;  
}

void mutgene_set_gene_name(MutGene *g, GtStr *gn) {
  gt_assert(g);
  g->gene_name = gn;  
}

GtStr* mutgene_get_type(MutGene *g) {
  gt_assert(g);
  return g->type;  
}

void mutgene_set_type(MutGene *g, GtStr *t) {
  gt_assert(g);
  g->type = t;  
}

unsigned long mutgene_get_rng_start(MutGene *g) {
  gt_assert(g);
  return g->rng_start;  
}

void mutgene_set_rng_start(MutGene *g, unsigned long rs) {
  gt_assert(g);
  g->rng_start = rs;  
}

unsigned long mutgene_get_rng_end(MutGene *g) {
  gt_assert(g);
  return g->rng_end;  
}


void mutgene_set_rng_end(MutGene *g, unsigned long re) {
  gt_assert(g);
  g->rng_end = re;  
}

unsigned long mutgene_get_phase(MutGene *g) {
  gt_assert(g);
  return g->phase;  
}

void mutgene_set_phase(MutGene *g, unsigned long p) {
  gt_assert(g);
  g->phase = p;  
}

unsigned long mutgene_get_child_size(MutGene *g) {
  gt_assert(g);
  return g->child_size;  
}

void mutgene_set_child_size(MutGene *g, unsigned long s) {
  gt_assert(g);
  g->phase = s;  
}

GtArray* mutgene_get_children_array(MutGene *g) {
  gt_assert(g);
  return g->children;  
}

void mutgene_add_child(MutGene *g, MutGene *child) {
  //g->children = gt_array_new(sizeof (child));
  gt_array_add_elem(g->children, child, sizeof(MutGene));
  g->child_size++;
  //~ printf("%lu",g->child_size);
}

void mutgene_add_content(MutGene *g, GtStr *i, GtStr *gn, GtStr *t, unsigned long rs, unsigned long re, unsigned long p) {
  mutgene_set_id(g,i);
  mutgene_set_gene_name(g,gn);
  mutgene_set_type(g,t);
  mutgene_set_rng_start(g,rs);
  mutgene_set_rng_end(g,re);
  mutgene_set_phase(g,p);  
}


void mutgene_reset(MutGene *g) {
  gt_assert(g);
}

void mutgene_delete(MutGene *g) {
  //~ unsigned long i = 0;
  gt_free(g->id);
  gt_free(g->gene_name);
  gt_free(g->type);  
  //~ for(i=0;i<gt_array_size(mutgene_get_children_array(g));i++){ 
    //~ mutgene_delete(gt_array_get(mutgene_get_children_array(g),i));
  //~ }
  gt_free(g->children);
  gt_free(g);
}