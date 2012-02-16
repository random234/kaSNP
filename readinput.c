#include "genometools.h"
#include "readinput.h"
#include "gff3vis.h"
#include "core/str_api.h"
#include "core/io_api.h"
#include "core/tokenizer_api.h"
//#define DEBUG

int read_input(const char * gff_file, const char *vcf_file_input, const char *encseq_file_input, const char *out_file_input, unsigned long splice_site_range) {
  GtNodeStream *in_stream, *visitor_stream;
  GtNodeVisitor *gff3vis;
  gt_lib_init(); /* we use the libgenometools.so */
  GtIO *vcf_io;
  GtTokenizer *vcf_token;
  GtStr *encseq_file;
  GtError *err; 
  
  err = gt_error_new();
    
  encseq_file = gt_str_new_cstr(encseq_file_input);
  gt_assert(encseq_file);
  
  vcf_io = gt_io_new(vcf_file_input, "r");
  gt_assert(vcf_io);
  
  vcf_token = gt_tokenizer_new(vcf_io);
  gt_assert(vcf_token);
  
  gt_tokenizer_skip_comment_lines(vcf_token); /* we skip the comment lines of the vcf file */
    
  
  gff3vis = gt_gff3_feat_vis_new(vcf_token, encseq_file, out_file_input, splice_site_range);
  in_stream = gt_gff3_in_stream_new_sorted(gff_file);
  gt_assert(in_stream);
  
  visitor_stream = gt_visitor_stream_new(in_stream, gff3vis);
  gt_assert(visitor_stream);

  gt_node_stream_pull(visitor_stream, err);  
  
  
  gt_node_stream_delete(visitor_stream);
  gt_node_stream_delete(in_stream);
  gt_tokenizer_delete(vcf_token);
    
  
  if (gt_lib_clean())
    return GT_EXIT_PROGRAMMING_ERROR; /* programmer error */
  return EXIT_SUCCESS;    
}

