#ifndef VCFOUTPUT_H_
#define VCFOUTPUT_H_

typedef struct VcfOutput VcfOutput;

VcfOutput* vcfoutput_new(void);

void vcfoutput_init(VcfOutput *, GtStr *);

void vcfoutput_write(VcfOutput *, ResultSet *);
 

#endif