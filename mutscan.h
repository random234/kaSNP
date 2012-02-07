#ifndef MUTSCAN_H
#define MUTSCAN_H
#include "mutgene.h"

typedef struct MutScan MutScan;


MutScan* mutscan_new(void);

GtStrArray* mutscan_get_vcf_array(MutScan *);

void mutscan_set_vcf_array(MutScan *, GtStrArray *);

MutGene* mutscan_get_mut_gene(MutScan *);

void mutscan_set_mut_gene(MutScan *, MutGene *);

unsigned long mutscan_get_splice_site_interval(MutScan *);

void mutscan_set_splice_site_interval(MutScan *, unsigned long);

unsigned long mutscan_init(MutScan *, GtStrArray *, GtFeatureNode *);

unsigned long mutscan_start_scan(MutScan *);

/* frame shift detection */
GtStrArray* mutscan_frms(MutScan *);

/* missense nonsense mutation detection */
GtStrArray* mutscan_miss(MutScan *mut);

/* splice site detection */
GtStrArray* mutscan_splice(MutScan *mut);

/* 5' UTR 2000bp Upstream detection 3' 500bp downstream detection */

/* intron mutation detection */
GtStrArray* mutscan_intron(MutScan *mut);



void mutscan_reset(MutScan *mut);

void mutscan_delete(MutScan*);

#endif