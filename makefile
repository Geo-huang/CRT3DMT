#-----------------------------------Part one-----------------------------------------
CXX                           := icpc
#CXXFLAGS                      := -std=c++11 -O3 -w -g -Wfatal-errors -heap -arrays -diag-disable=10441
CXXFLAGS                      := -std=c++11 -O3 -w -g -Wfatal-errors -heap -arrays -diag-disable=10441 -qopenmp 
EIGEN_INCLUDE                 := -I$(PWD)/contrib/eigen3.3.9
PARDISO_LIB                   := -mkl # automatically linking all used paralleled MKL libraries 


# For compiling
INCLUDES                      := -I$(PWD)/ -I$(PWD)/src/ -I$(PWD)/src/MT3D/
# For head files checking
INCLUDE_H                     := $(wildcard src/*.h) $(wildcard src/MT3D/*.h) 
# For compiling and checking
SRCS                          := $(wildcard *.cpp) $(wildcard src/*.cpp) $(wildcard src/MT3D/*.cpp)
# The directory of objects
OBJDIR                        := ./src/OBJS
OBJS                          := $(patsubst %.cpp, $(OBJDIR)/%.o, $(SRCS))


#-----------------------------------Part two-----------------------------------------

# Here is the compile step
$(OBJDIR)/%.o: %.cpp 
	@echo "CRT3DMT is compiling C++ "$<"..."
	@$(CXX) $(CXXFLAGS) $(INCLUDES) $(EIGEN_INCLUDE) -c $< -o $@


# Only do compilation
all: $(OBJS)


#-----------------------------------Part three---------------------------------------

# Here is the link step	
CRT3DMT: $(OBJS)  
	 @$(CXX) $(CXXFLAGS) $(PARDISO_LIB) $(OBJS) -o CRT3DMT 


#-----------------------------------Part four----------------------------------------

# Type "make clean" to get rid of all objects and executable files
clean: 
	@rm -rf $(OBJS) CRT3DMT 


#-----------------------------------Part five-----------------------------------------

# Check the compile and link options.	
echo:	
	@echo "C++ Compiler:                         \n$(CXX)\n"
	@echo "CXXFLAGS:                             \n$(CXXFLAGS)\n"
        # Check head files 
	@echo "Header Files in src folder:           \n$(INCLUDE_H)\n"
        # Check source files
	@echo "Source Files in src folder:           \n$(SRCS)\n"
        # Check object files
	@echo "Object Files in src folder:           \n$(OBJS)\n"
	
	
