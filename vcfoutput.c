#include <string.h>
#include "genometools.h"
#include "core/file_api.h"
#include "core/str_array_api.h"
#include "core/str_api.h"
#include "resultset.h"
#include "vcfoutput.h"

struct VcfOutput{
  GtFile *file;
};

VcfOutput* vcfoutput_new(void) {
  VcfOutput *v = gt_malloc(sizeof(*v));
  return v;
}

void vcfoutput_init(VcfOutput *v, const char *of) {
  gt_assert(v);
  gt_assert(of);
  GtError *err = gt_error_new();
  v->file = gt_file_new(of,"w", err);
}

void vcfoutput_write(VcfOutput *v, ResultSet *r) {
  gt_assert(v);
  gt_assert(r); 
  GtStr *temp = gt_str_new();
  unsigned long i =0;
  unsigned long vcf_size = 0;
  
  vcf_size = gt_str_array_size(resultset_get_vcf_array(r));
  for(i=0;i<gt_str_array_size(resultset_get_vcf_array(r));i++) {
    gt_str_set(temp, gt_str_array_get(resultset_get_vcf_array(r),i));
    
    if(i == vcf_size-1) {
      gt_str_append_cstr(temp,";");
      
      if(resultset_get_exon(r) != 0) {
        gt_str_append_cstr(temp,"EX;");
      }
      if(resultset_get_frms(r) != 0) {
        gt_str_append_cstr(temp,"NSF;");
      }
      if(resultset_get_miss(r) != 0) {
        gt_str_append_cstr(temp,"NSM;");
      }
      if(resultset_get_nons(r) != 0) {
        gt_str_append_cstr(temp,"NSN;");
      }      
      if(resultset_get_threeprime(r) != 0) {
        gt_str_append_cstr(temp,"ASS;");
      }
      if(resultset_get_fiveprime(r) != 0) {
        gt_str_append_cstr(temp,"DSS;");
      }
    }
    gt_str_append_cstr(temp,"\t");
    gt_file_xwrite(v->file,gt_str_get(temp),gt_str_length(temp));    
    gt_str_reset(temp);
  } 
  gt_str_delete(temp);
}

void vcfoutput_delete(VcfOutput *v){
  gt_file_delete(v->file);
  
}