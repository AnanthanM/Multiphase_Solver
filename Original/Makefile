NAME 	:= fvmvof
OBDIR	:= OBJ/
RUNDIR	:= RUN

OBJECTS =  $(OBDIR)advection.o   $(OBDIR)fvmvof.o $(OBDIR)diffusion.o $(OBDIR)pressure.o $(OBDIR)boundary.o

HEADER 	:=	fvmvof.h\
		boundary.h

.DEFAULT:

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm 

CC	:=	gcc
FLAGS	:=	  $(INCPATH) -Wall -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=1 
#-fopenmp
#FLAGS	:=	-g -pg $(INCPATH) -Wall -fopenmp -O3
#LFLAGS	:=	$(LIBPATH) $(LIBS) -lm -lgsl -lblas
LFLAGS	:=	$(LIBPATH) $(LIBS) 

####################################################################

$(NAME): $(OBJECTS) $(RUNDIR)
	$(CC) -o $(NAME) $(OBJECTS) $(FLAGS) $(LFLAGS)
	@mv $(NAME) $(RUNDIR) 
#	@cp iwrite $(RUNDIR)

$(RUNDIR): 
	@test -d $(RUNDIR) || mkdir $(RUNDIR)

$(OBDIR)fvmvof.o: fvmvof.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)fvmvof.o fvmvof.c

$(OBDIR)advection.o: advection.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)advection.o advection.c

$(OBDIR)diffusion.o: diffusion.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)diffusion.o diffusion.c

$(OBDIR)pressure.o: pressure.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)pressure.o pressure.c

$(OBDIR)boundary.o: boundary.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)boundary.o boundary.c


copy:
	cp -f clean $(RUNDIR)/.	
clean: 		
	rm -f *~
	rm $(OBDIR)*.o
