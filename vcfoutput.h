#ifndef VCFOUTPUT_H_
#define VCFOUTPUT_H_

typedef struct VcfOutput VcfOutput;

VcfOutput* vcfoutput_new(void);

void vcfoutput_init(VcfOutput *, const char *);

void vcfoutput_write(VcfOutput *, ResultSet *);
 
void vcfoutput_delete(VcfOutput *);

#endif