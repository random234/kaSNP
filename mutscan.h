#ifndef MUTSCAN_H
#define MUTSCAN_H

typedef struct MutScan MutScan;


MutScan* mutscan_new(void);

GtStrArray* mutscan_frms(MutScan *, GtStrArray *, GtFeatureNode *);

GtStrArray* mutscan_miss(MutScan *mut);

GtStrArray* mutscan_nons(MutScan *mut);

void mutscan_reset(MutScan *mut);

void mutscan_delete(MutScan*);

#endif