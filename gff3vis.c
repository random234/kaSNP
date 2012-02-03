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
#include "gff3vis.h"
#include "mutscan.h"


/*
#include "vcfentry.h"
#include "mutscan.h"
#include "gff3data.h"
#include "gff3entry.h"
*/


struct GtGff3Vis {
  const GtNodeVisitor parent_instance;
  GtTokenizer *vcf_token;
  GtStr *fas_file;
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
  GtFeatureNodeIterator *fni, *fni_childs;
  GtFeatureNode *node, *child;
  GtGff3Vis *v;
  int had_err = 0;
  GtSplitter *vcf_split;
  GtStr *line, *tokenline;
  GtStrArray *vcf_arr;
  GT_UNUSED MutScan *mutscan;  
  unsigned long pos = 0;
  unsigned long i;
  gt_error_check(err);
  v = gt_gff3_vis_cast(nv);

  /* no more lines in VCF */
  if (!gt_tokenizer_has_token(v->vcf_token))
    return 0;
    
  fni = gt_feature_node_iterator_new(fn);
  vcf_split = gt_splitter_new();
  line = gt_str_new();
  tokenline = gt_tokenizer_get_token(v->vcf_token);
  gt_str_append_str(line, tokenline);
  gt_str_delete(tokenline);
  gt_splitter_split(vcf_split, gt_str_get(line),
                    gt_str_length(line), '\t');
  GT_SNPANNO_CHECK_FILE_FORMAT;

  while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
    GtRange rng = gt_genome_node_get_range((GtGenomeNode*) node);    
    while(!had_err && pos <= rng.end) {
      /* printf("beg splitter size: %lu\n", gt_splitter_size(vcf_split)); */
      pos = 0; pos = strtol(gt_splitter_get_token(vcf_split,1),NULL,0);
      
      if(strcmp(gt_feature_node_get_type(node), "gene") == 0 && rng.start <= pos && rng.end >= pos && 
        (strcmp(gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)),gt_splitter_get_token(vcf_split,0)) == 0)){

        printf("SNP in gene found at Pos: %lu\t",pos);
        printf("rng.end: %lu \n",rng.end);
          
        printf("type: %s, %lu-%lu, seqid %s\n\n",
                gt_feature_node_get_type(node),
                rng.start,
                rng.end,
                gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) node)));
      
        fni_childs =  gt_feature_node_iterator_new_direct(node);      
        while ((child = gt_feature_node_iterator_next(fni_childs))) {
          //~ GtRange rng_child = gt_genome_node_get_range((GtGenomeNode*) child);
          //~ printf("CHILD: type: %s, %lu-%lu, seqid %s\n",
                //~ gt_feature_node_get_type(child),
                //~ rng_child.start,
                //~ rng_child.end,
                //~ gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) child)));        
        }        
        gt_feature_node_iterator_delete(fni_childs);

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
        /* put VCF entries into gt_str_array object and pass them to MutScan object */
        /* VCF entry: CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO */
        vcf_arr = gt_str_array_new();
        for(i=0;i<gt_splitter_size(vcf_split);i++) {
          gt_str_array_add_cstr(vcf_arr,gt_splitter_get_token(vcf_split, i));
        }
        
        
        
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
  gt_splitter_delete(vcf_split);
  return had_err;
}





void gt_gff3_vis_free(GtNodeVisitor *nv)
{
  GtGff3Vis *v;
  if (!nv) return;
  v = gt_gff3_vis_cast(nv);
  //gt_str_array_delete(v->sa);
  gt_str_delete(v->fas_file);
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

GtNodeVisitor* gt_gff3_feat_vis_new(GtTokenizer *vcf_token, GtStr *fas_file)
{
  GtGff3Vis *mfv;
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(gt_gff3_vis_class());
  mfv = gt_gff3_vis_cast(nv);
  //mfv->sa = gt_str_array_new();
  mfv->vcf_token = vcf_token;
  mfv->fas_file = fas_file;  
  return nv;
}

//~ static int gt_gff3_vis_feature_node(GtNodeVisitor *nv,
                                           //~ GtFeatureNode *fn,
                                            //~ GtError *err)
//~ {
  //~ GtFeatureNodeIterator *fni, *fni_childs;
  //~ GtFeatureNode *node, *child;
  //~ GtGff3Vis *v;
  //~ GtSplitter *vcf_split;  
  //~ GtRange rng;
  //~ gt_error_check(err);
  //~ v = gt_gff3_vis_cast(nv);
  //~ fni = gt_feature_node_iterator_new(fn);
  //~ unsigned long len = 0;
  //~ unsigned long vcf_pos = 0;
    

    
  //~ vcf_split = gt_splitter_new();
  //~ len = gt_str_length(gt_tokenizer_get_token(v->vcf_token)); // gt_str_length
  //~ gt_splitter_split(vcf_split, gt_str_get(gt_tokenizer_get_token(v->vcf_token)),len, '\t');  
  //~ vcf_pos = 0; vcf_pos = strtol(gt_splitter_get_token(vcf_split,1),NULL,0);
  //~ rng.start = 0; rng.end = 1;
  //~ while((node = gt_feature_node_iterator_next(fni))) {
    
    //~ // wenn gen gefunden verarbeite es ansonsten hole naechsten node 
    //~ if(strcmp(gt_feature_node_get_type(node), "gene") == 0) {
      //~ rng = gt_genome_node_get_range((GtGenomeNode*) node);
      //~ printf("found gene\n");
      
      //~ // solange aktuelle vcf_pos kleiner als start des aktuell betrachteten Gens hole naechsten VCF Eintrag
      //~ while(vcf_pos<rng.start) {
        //~ printf("increasing VCF position because rng.start<vcf_pos POS: %lu rng.start: %lu\n", vcf_pos, rng.start);
        //~ gt_tokenizer_next_token(v->vcf_token);
        //~ gt_splitter_reset(vcf_split);
        //~ len=0;  len = gt_str_length(gt_tokenizer_get_token(v->vcf_token));
        //~ gt_splitter_split(vcf_split, gt_str_get(gt_tokenizer_get_token(v->vcf_token)),len, '\t');
        //~ vcf_pos = 0; vcf_pos = strtol(gt_splitter_get_token(vcf_split,1),NULL,0);
      //~ }
      
      
      
      //~ // wenn gen gefunden und vcf_pos innerhalb des gens hole alle Kindelemente des zugehoerigen Gens
      //~ while(vcf_pos<=rng.end && vcf_pos >= rng.start){   
        //~ fni_childs =  gt_feature_node_iterator_new_direct(node);      
        //~ while ((child = gt_feature_node_iterator_next(fni_childs))) {
          //~ printf("getting child elements for VCF entry POS: %lu\n",vcf_pos);
          
          
          
          
        //~ }
        //~ VcfEntry *vcf_entry = vcf_entry_new(vcf_split);
        //~ Gff3Data *gff3data;
        //~ check_frame(gff3data, vcf_entry);
        
        //~ gt_feature_node_iterator_delete(fni_childs);
        
        //~ gt_tokenizer_next_token(v->vcf_token);
        //~ gt_splitter_reset(vcf_split);
        //~ len=0;  len = gt_str_length(gt_tokenizer_get_token(v->vcf_token));
        //~ gt_splitter_split(vcf_split, gt_str_get(gt_tokenizer_get_token(v->vcf_token)),len, '\t');
        //~ vcf_pos = 0; vcf_pos = strtol(gt_splitter_get_token(vcf_split,1),NULL,0);
      //~ }

      //~ while(vcf_pos>rng.end) {
        //~ printf("increasing GFF3 position because vcf_pos>rng.end POS: %lu rng.end: %lu\n", vcf_pos, rng.end);
        //~ node = gt_feature_node_iterator_next(fni);
        //~ rng = gt_genome_node_get_range((GtGenomeNode*) node);
      //~ }
    //~ }
  //~ }
  //~ gt_feature_node_iterator_delete(fni);
  //~ return 0;
//~ }

