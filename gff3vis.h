#ifndef GFF3VIS_H
#define GFF3VIS_H
#include "core/tokenizer_api.h"


/* Implements the <GtNodeVisitor> interface. */
//typedef struct GtMyFeatureVisitor GtMyFeatureVisitor;
typedef struct GtGff3Vis GtGff3Vis;


#include "extended/node_visitor_api.h"

//const GtNodeVisitorClass* gt_my_feature_visitor_class(void);
//GtNodeVisitor*            gt_my_feature_visitor_new(void);

const GtNodeVisitorClass* gt_gff3_vis_class(void);
GtNodeVisitor* gt_gff3_feat_vis_new(GtTokenizer *, GtStr *, const char *, unsigned long);
unsigned long get_description_file_number(GtStrArray *, GtStr *);

#endif
