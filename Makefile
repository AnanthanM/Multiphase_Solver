NAME 	:= EXECUTABLE
OBDIR	:= OBJ/
RUNDIR	:= RUN

OBJECTS =  $(OBDIR)Field_Allocation.o  $(OBDIR)VTK_Output.o $(OBDIR)Set_Ghost_Cells.o  $(OBDIR)FVM.o 

HEADER 	:=	fvm.h

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm 

CC	:=	gcc
FLAGS	:=	  $(INCPATH) 
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

$(OBDIR)Field_Allocation.o: Field_Allocation.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Field_Allocation.o Field_Allocation.c

clean: 		
	rm -f *~
	rm $(OBDIR)*.o

Delete_Results: 		
	rm -f *~
	rm $(RUNDIR)/*.vtk
	rm $(RUNDIR)/$(NAME)
