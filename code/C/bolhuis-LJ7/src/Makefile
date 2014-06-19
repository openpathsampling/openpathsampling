SRCS=	ran.c main.c  mcsteps.c md.c construct.c draws.c 

OBJS=	ran.c main.c  mcsteps.c md.c construct.c draws.c


HFILES = 

#

CF =   cc
FFLAGS =   
LIBS =   -lm -framework OpenGL -framework Glut -framework Cocoa
OUT =	path.run



opt:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -o $(OUT) $(SRCS:.c=.o) $(LIBS)

prof:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -pg -o $(OUT) $(SRCS) $(LIBS)

opt3:	$(SRCS)
	$(CF)  -O3 -o $(OUT) $(SRCS) $(LIBS)

deb:	$(SRCS)
	$(CF) -g -O2 -o rods.deb $(SRCS) $(LIBS)

$(SRCS):	$(HFILES)
$(OBJS):	$(HFILES)	







