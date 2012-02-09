#ifndef MUTSCAN_H
#define MUTSCAN_H
#include "mutgene.h"
#include "resultset.h"

typedef struct MutScan MutScan;


MutScan* mutscan_new(void);

GtStrArray* mutscan_get_vcf_array(MutScan *);

void mutscan_set_vcf_array(MutScan *, GtStrArray *);

MutGene* mutscan_get_mut_gene(MutScan *);

void mutscan_set_mut_gene(MutScan *, MutGene *);

unsigned long mutscan_get_splice_site_interval(MutScan *);

void mutscan_set_splice_site_interval(MutScan *, unsigned long);

unsigned long mutscan_init(MutScan *, GtStrArray *, GtFeatureNode *,GtEncseq *,GtEncseqReader *);

ResultSet* mutscan_start_scan(MutScan *);

/* frame shift detection */
unsigned long mutscan_frame(MutScan *, ResultSet*);

/* missense nonsense mutation detection */
unsigned long mutscan_miss(MutScan *mut, ResultSet*);

/* splice site detection */
unsigned long mutscan_splice(MutScan *mut, ResultSet*);

/* 5' UTR 2000bp Upstream detection 3' 500bp downstream detection */
unsigned long mutscan_utr(MutScan *mut, ResultSet*);

/* intron mutation detection */
unsigned long mutscan_exon(MutScan *mut, ResultSet*);


void mutscan_reset(MutScan *mut);

void mutscan_delete(MutScan*);

#endif