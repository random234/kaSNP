#ifndef MUTSCAN_H
#define MUTSCAN_H

typedef struct MutScan MutScan;


MutScan* mutscan_new(void);

unsigned long mutscan_init(MutScan *, GtStrArray *, GtFeatureNode *);

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