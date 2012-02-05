#ifndef MUTGENE_H
#define MUTGENE_H

typedef struct MutGene MutGene;

MutGene* mutgene_new(void);

void mutgene_add_content(MutGene *,GtStr *, GtRange *,unsigned long);

void mutgene_add_child(MutGene *, MutGene *);


void mutgene_reset(MutGene*);

void mutgene_delete(MutGene*);

#endif