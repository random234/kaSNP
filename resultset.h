#ifndef RESULTSET_H_
#define RESULTSET_H_

typedef struct ResultSet ResultSet;

ResultSet* resultset_new(void);

void resultset_set_vcf_array(ResultSet *, GtStrArray *);

GtStrArray * resultset_get_vcf_array(ResultSet *);

void resultset_set_var_pos(ResultSet *, unsigned long );

unsigned long resultset_get_var_pos(ResultSet *);

void resultset_set_id(ResultSet *, GtStr *);

GtStr * resultset_get_id(ResultSet *);

void resultset_set_gene_name(ResultSet *, GtStr *);

GtStr * resultset_get_gene_name(ResultSet *);

void resultset_add_mrna_id(ResultSet *, GtStr *);

GtStrArray * resultset_get_mrna_ids(ResultSet *);

unsigned long resultset_check_mrna_ids(ResultSet *, GtStr *);

void resultset_set_frms(ResultSet *, unsigned long);
  
unsigned long resultset_get_frms(ResultSet *);

void resultset_add_dna_seq(ResultSet *, GtStr *);

GtStrArray * resultset_get_dna_seqs(ResultSet *);
  
void resultset_set_miss(ResultSet *, unsigned long);

unsigned long resultset_get_miss(ResultSet *);

void resultset_set_nons(ResultSet *, unsigned long);

unsigned long resultset_get_nons(ResultSet *);

#endif