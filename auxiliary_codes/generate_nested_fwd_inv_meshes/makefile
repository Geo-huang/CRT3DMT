#Use value and fuction
SHELL        := /bin/bash
CXX          := icpc
CXXFLAGS     := -std=c++11 -O3 -g 

INCLUDES     := $(wildcard /*.h)
SRCS         := $(wildcard *.cpp)
OBJS         := $(patsubst %.cpp,%.o,$(SRCS))



#How to compile C++
%.o: %.cpp
	@echo "MakeNestedMesh is compling: "$<"..."
	@$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< -o $@


#How to link C++
MakeNestedMesh: $(OBJS)
	@$(CXX) $(CXXFLAGS) -o MakeNestedMesh $(OBJS) 
	@chmod a+x MakeNestedMesh


#What to clean
clean:
	@rm -rf *.o MakeNestedMesh


echo:	
	@echo -e "Header Files:\n$(INCLUDES)\n"
	@echo -e "Source Files:\n$(SRCS)\n"
	@echo -e "Object Files:\n$(OBJS)\n"
	@echo -e "C++ Compiler:\n$(CXX)\n"
	@echo -e "CXXFLAGS:    \n$(CXXFLAGS)\n"



