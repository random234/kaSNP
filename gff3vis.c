#include <string.h>
#include "genometools.h"
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/ma_api.h"
#include "core/str_array_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "core/splitter_api.h"
#include "core/str_array_api.h"
#include "core/encseq_api.h"
#include "core/readmode_api.h"
#include "gff3vis.h"
#include "mutscan.h"
#include "vcfoutput.h"



struct GtGff3Vis {
  const GtNodeVisitor parent_instance;
  GtTokenizer *vcf_token;
  GtEncseq *encseq;
  VcfOutput *vcf_out;
  GtEncseqLoader *encseq_load;
    
  unsigned long splice_site_range;
  //GtStrArray *sa;
};

#define gt_gff3_vis_cast(GV)\
        gt_node_visitor_cast(gt_gff3_vis_class(), GV)

#define GT_SNPANNO_CHECK_FILE_FORMAT \
  if (!had_err && gt_splitter_size(vcf_split) != 8) { \
    had_err = -1; \
    gt_error_set(err, "found %lu columns in file %s, line %lu instead of 8!", \
                      gt_splitter_size(vcf_split), \
                      gt_tokenizer_get_filename(v->vcf_token), \
                      gt_tokenizer_get_line_number(v->vcf_token)); \
  } \

static int gt_gff3_vis_feature_node(GtNodeVisitor *nv,
                                           GtFeatureNode *fn,
                                            GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  GtGff3Vis *v;
  GtEncseqReader *encseq_read;
  int had_err = 0;
  GtSplitter *vcf_split;
  GtSplitter *desc_split;
  GtStr *line, *tokenline;
  GtStrArray *vcf_arr, *desc_arr;
  MutScan *mut;  
  unsigned long pos = 0;
  unsigned long i;
  unsigned long desclen, chrom_file;
  GtStr *temp;
  gt_error_check(err);
  v = gt_gff3_vis_cast(nv);
  
  encseq_read = gt_encseq_create_reader_with_readmode(v->encseq,0, 0);
  
  /* no more lines in VCF */
  if (!gt_tokenizer_has_token(v->vcf_token))
    return 0;
    
  fni = gt_feature_node_iterator_new(fn);
  vcf_arr = gt_str_array_new();
  desc_arr = gt_str_array_new();
  temp = gt_str_new();
  mut = mutscan_new();
  vcf_split = gt_splitter_new();
  desc_split = gt_splitter_new();
  line = gt_str_new();
  tokenline = gt_tokenizer_get_token(v->vcf_token);
  gt_str_append_str(line, tokenline);
  gt_str_delete(tokenline);
  gt_splitter_split(vcf_split, gt_str_get(line),
                    gt_str_length(line), '\t');
  GT_SNPANNO_CHECK_FILE_FORMAT;
  
  /* we split the description file in to access the first field during mutation scanning */
  for(i=0;i< gt_str_array_size(gt_encseq_filenames(v->encseq));i++){    
    temp = gt_str_new_cstr(gt_encseq_description(v->encseq,&desclen,i));       
    gt_splitter_split(desc_split, gt_str_get(temp),strlen(gt_str_get(temp)), ' ');  
    gt_str_array_add(desc_arr,gt_str_new_cstr(gt_splitter_get_token(desc_split,0)));
    //~ printf("desc: %s \n", gt_str_array_get(desc_arr,i));
    gt_splitter_reset(desc_split);
    gt_str_reset(temp);
  }
  
  chrom_file = get_description_file_number(desc_arr,gt_genome_node_get_seqid((GtGenomeNode*) fn));
  
  
  while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
    GtRange rng = gt_genome_node_get_range((GtGenomeNode*) node);      
    
    /* check if seqid is equal to current file_mapping for the encseq access */
    if(strcmp(gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)), gt_str_array_get(desc_arr,chrom_file)) != 0) {
      chrom_file = get_description_file_number(desc_arr,gt_genome_node_get_seqid((GtGenomeNode*) fn));
    }
    
    while(!had_err && pos <= rng.end) {
      /* printf("beg splitter size: %lu\n", gt_splitter_size(vcf_split)); */
      pos = 0; pos = strtol(gt_splitter_get_token(vcf_split,1),NULL,0);
      
      if(strcmp(gt_feature_node_get_type(node), "gene") == 0 && rng.start <= pos && rng.end >= pos && 
        (strcmp(gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)),gt_splitter_get_token(vcf_split,0)) == 0)){

        //~ printf("SNP in gene found at Pos: %lu\t",pos);
        //~ printf("rng.end: %lu \n",rng.end);
          
        //~ printf("type: %s, %lu-%lu, seqid %s\n\n",
                //~ gt_feature_node_get_type(node),
                //~ rng.start,
                //~ rng.end,
                //~ gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)));
      
        gt_tokenizer_next_token(v->vcf_token);
        if (!gt_tokenizer_has_token(v->vcf_token))
          break;
        gt_splitter_reset(vcf_split);
        gt_str_reset(line);
        tokenline = gt_tokenizer_get_token(v->vcf_token);
        gt_str_append_str(line, tokenline);
        gt_str_delete(tokenline);
        gt_splitter_split(vcf_split, gt_str_get(line),
                          gt_str_length(line), '\t');
        GT_SNPANNO_CHECK_FILE_FORMAT;
        /* put VCF entries into gt_str_array object */
        /* VCF entry: CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO */
        gt_str_array_reset(vcf_arr);
        vcf_arr = gt_str_array_new();
        for(i=0;i<gt_splitter_size(vcf_split);i++) {
          gt_str_array_add_cstr(vcf_arr,gt_splitter_get_token(vcf_split, i));
        }
        /* pass VCF entries and FeatureNode containing child nodes to MutScan object */
        //~ if(!mutscan_init(mut, vcf_arr, node))
          //~ break;
        
        //~ printf("%lu",v->splice_site_range);
        mutscan_init(mut, vcf_arr, chrom_file, node, v->encseq, encseq_read);
        vcfoutput_write(v->vcf_out,mutscan_start_scan(mut));
        
        mutscan_reset(mut);
        
      } else {        
        //~ printf("SNP not in gene at Pos: %lu\t",pos);
        //~ printf("rng.end: %lu\n \n",rng.end);
        
        gt_tokenizer_next_token(v->vcf_token);
        if (!gt_tokenizer_has_token(v->vcf_token))
          break;
        gt_splitter_reset(vcf_split);
        gt_str_reset(line);
        tokenline = gt_tokenizer_get_token(v->vcf_token);
        gt_str_append_str(line, tokenline);
        gt_str_delete(tokenline);
        gt_splitter_split(vcf_split, gt_str_get(line),
                          gt_str_length(line), '\t');
        GT_SNPANNO_CHECK_FILE_FORMAT;
      }
    }
  }
  gt_feature_node_iterator_delete(fni);
  gt_str_delete(line);
  gt_str_delete(temp);
  gt_str_array_delete(vcf_arr);
  gt_str_array_delete(desc_arr);
  gt_splitter_delete(vcf_split);
  gt_splitter_delete(desc_split);
  mutscan_delete(mut);  
  return had_err; 
}

unsigned long get_description_file_number(GtStrArray *desc, GtStr *seqid) {
  unsigned long i = 0;  
  for(i=0;i<gt_str_array_size(desc);i++) {
    if(strcmp(gt_str_array_get(desc,i), gt_str_get(seqid)) == 0) {
      //~ printf("BALRG: %lu \n", i);
      return i;
    } 
  }
  return 0;
}



void gt_gff3_vis_free(GtNodeVisitor *nv)
{
  GtGff3Vis *v;
  if (!nv) return;
  v = gt_gff3_vis_cast(nv);
  vcfoutput_delete(v->vcf_out);
  //gt_str_array_delete(v->sa);
  //~ gt_str_delete(v->fas_file);
}

const GtNodeVisitorClass* gt_gff3_vis_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtGff3Vis),
                                    gt_gff3_vis_free,
                                    NULL,
                                    gt_gff3_vis_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_gff3_feat_vis_new(GtTokenizer *vcf_token, GtStr *encseq_file, const char *out_file, unsigned long splice_site_range)
{
  GtError *err;
  err = gt_error_new();
  GtGff3Vis *mfv;
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(gt_gff3_vis_class());
  mfv = gt_gff3_vis_cast(nv);
  mfv->vcf_out = vcfoutput_new();
  vcfoutput_init(mfv->vcf_out, out_file);
  mfv->vcf_token = vcf_token;
  mfv->encseq_load = gt_encseq_loader_new();
  mfv->encseq = gt_encseq_loader_load(mfv->encseq_load,gt_str_get(encseq_file), err);
  mfv->splice_site_range = splice_site_range;
  return nv;
}

