#include "genometools.h"
#include "mutscan.h"

struct MutScan {
  char lala;
  
  
};

MutScan* mutscan_new(void) {
  MutScan *mut = gt_malloc(sizeof (MutScan));
  
  return mut;
}

void mutscan_delete(MutScan *mut) {
  gt_free(mut);  
}