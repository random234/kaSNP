#ifndef MUTGENE_H
#define MUTGENE_H

typedef struct MutGene MutGene;

MutGene* mutgene_new(void);

GtStr* mutgene_get_id(MutGene *);

void mutgene_set_id(MutGene *, GtStr *);

GtStr* mutgene_get_gene_name(MutGene *);

void mutgene_set_gene_name(MutGene *, GtStr *);

GtStr* mutgene_get_type(MutGene *);

void mutgene_set_type(MutGene *, GtStr *);

unsigned long mutgene_get_rng_start(MutGene*);

void mutgene_set_rng_start(MutGene*, unsigned long);

unsigned long mutgene_get_rng_end(MutGene*);

void mutgene_set_rng_end(MutGene*, unsigned long);

unsigned long mutgene_get_phase(MutGene*);

void mutgene_set_phase(MutGene*, unsigned long);

unsigned long mutgene_get_child_size(MutGene*);

void mutgene_set_child_size(MutGene*, unsigned long);

GtArray* mutgene_get_children_array(MutGene *);

void mutgene_add_child(MutGene*, MutGene*);

void mutgene_add_content(MutGene*,GtStr*, GtStr*, GtStr*, unsigned long, unsigned long, unsigned long);

void mutgene_reset(MutGene*);

void mutgene_delete(MutGene*);

#endif