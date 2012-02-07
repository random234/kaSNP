# set prefix for libgenometools
prefix ?= /usr/local

CC:=gcc
# -Werror -Wall -Wunused-parameter
# -Os
GT_CFLAGS:=-Wall -Werror -Wunused-parameter -g -Os -pipe \
           -I$(prefix)/include/genometools \
           -I/usr/include/cairo \
           -I/usr/local/include/cairo \
           -I/sw/include/cairo \
           -I/opt/local/include/cairo
LD:=$(CC)
GT_LDFLAGS:=-L$(prefix)/lib \
            -L/usr/X11R6/lib \
            -L/usr/local/lib \
            -L/opt/local/lib

all: kasnp

kasnp: kasnp.o readinput.o gff3vis.o mutscan.o mutgene.o resultset.o
	$(LD) $(LDFLAGS) $(GT_LDFLAGS) $< -lm -lgenometools -lcairo -o $@ readinput.o gff3vis.o mutscan.o mutgene.o resultset.o

#%.o: %.c
#  gcc -Wall -g -c $<


# generic compilation rule which creates dependency file on the fly
%.o: %.c
	$(CC) -c $< -DWITHOUT_CAIRO -o $@ $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read dependencies
-include $(wildcard *.d)

.PHONY: clean
clean:
	rm -f *.[od]
