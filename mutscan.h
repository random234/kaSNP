#ifndef MUTSCAN_H
#define MUTSCAN_H

typedef struct MutScan MutScan;


MutScan* mutscan_new(void);
void mutscan_delete(MutScan*);

#endif