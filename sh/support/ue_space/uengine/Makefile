all:				ixp2000-ue-disas		\
				ixp2000-lib-userspace.o		\
				ixp2000-uengine.o		\
				ixp2400-msf.o

clean:
				rm -f ixp2000-ue-disas *.o

ixp2000-lib-userspace.o:	ixp2000-lib-userspace.c
				$(CC_CROSS) $(CFLAGS_ARM) -c -o ixp2000-lib-userspace.o ixp2000-lib-userspace.c

ixp2000-uengine.o:		ixp2000-uengine.c ixp2000-uengine.h
				$(CC_CROSS) $(CFLAGS_ARM) -c -o ixp2000-uengine.o ixp2000-uengine.c

ixp2000-ue-disas:		ixp2000-ue-disas.c ixp2000-ue-disas.h
				$(CC_CROSS) $(CFLAGS_ARM) -DTEST -o ixp2000-ue-disas ixp2000-ue-disas.c

ixp2400-msf.o:			ixp2400-msf.c ixp2400-msf.h
				$(CC_CROSS) $(CFLAGS_ARM) -c -o ixp2400-msf.o ixp2400-msf.c
