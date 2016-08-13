NAME 	:= EXECUTABLE
OBDIR	:= OBJ/
RUNDIR	:= 100

OBJECTS =  $(OBDIR)Helper_Functions.o $(OBDIR)VTK_Output.o $(OBDIR)Set_Ghost_Cells.o $(OBDIR)FVM.o $(OBDIR)VOF_Functions.o $(OBDIR)Advection.o $(OBDIR)Diffusion.o $(OBDIR)BiCGSTAB.o $(OBDIR)Pressure.o

HEADER 	:=	fvm.h

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm -lvofi

CC	:=	gcc 
FLAGS	:=	  $(INCPATH) -g
LFLAGS	:=	$(LIBPATH) $(LIBS) 

####################################################################

$(NAME): $(OBJECTS) $(RUNDIR)
	$(CC) -o $(NAME) $(OBJECTS) $(FLAGS) $(LFLAGS)
	@mv $(NAME) $(RUNDIR) 

$(RUNDIR): 
	@test -d $(RUNDIR) || mkdir $(RUNDIR)

$(OBDIR)FVM.o: FVM.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)FVM.o FVM.c

$(OBDIR)Set_Ghost_Cells.o: Set_Ghost_Cells.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Set_Ghost_Cells.o Set_Ghost_Cells.c

$(OBDIR)VTK_Output.o: VTK_Output.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)VTK_Output.o VTK_Output.c

$(OBDIR)VOF_Functions.o: VOF_Functions.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)VOF_Functions.o VOF_Functions.c

$(OBDIR)Helper_Functions.o: Helper_Functions.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Helper_Functions.o Helper_Functions.c

$(OBDIR)Advection.o: Advection.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Advection.o Advection.c

$(OBDIR)Diffusion.o: Diffusion.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Diffusion.o Diffusion.c

$(OBDIR)BiCGSTAB.o: BiCGSTAB.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)BiCGSTAB.o BiCGSTAB.c

$(OBDIR)Pressure.o: Pressure.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Pressure.o Pressure.c

clean: 		
	rm -f *~
	rm $(OBDIR)*.o

Delete_Results: 		
	rm -f *~
	rm $(RUNDIR)/*.vtk
	rm $(RUNDIR)/$(NAME)
