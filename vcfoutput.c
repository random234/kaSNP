#include <string.h>
#include "genometools.h"
#include "resultset.h"
#include "vcfoutput.h"

struct VcfOutput{
  GtStr *outfile;
  ResultSet *results;
};

VcfOutput* vcfoutput_new(void) {
  VcfOutput *v = gt_malloc(sizeof(*v));
  v->outfile = gt_str_new();
  return v;
}

void vcfoutput_init(VcfOutput *v, GtStr *of) {
  gt_assert(v);
  gt_assert(of);
  v->outfile = of;
}

void vcfoutput_write(VcfOutput *v, ResultSet *r) {
  gt_assert(v);
  gt_assert(r);
  v->results = r;
}