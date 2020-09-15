CC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11 

CHIMESFF_SRC=$(CC_LOC)/../../../../src/chimesFF.cpp
CHIMESFF_HDR=$(CC_LOC)/../../../../src/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)
		
WRAPPER_SRC=$(CC_LOC)/../src/wrapper-C.cpp
WRAPPER_HDR=$(CC_LOC)/../src/wrapper-C.h
		
wrapper-C.o : $(WRAPPER_SRC) $(WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(WRAPPER_SRC) $(CHIMESFF_SRC) -I $(CC_LOC)/../src/ -I $(CC_LOC)/../../../../src/
	
CC = gcc

test_cpp_and_c_wrapper.o : test_cpp_and_c_wrapper.c $(WRAPPER_SRC) $(WRAPPER_HDR)
	$(CC) -c test_cpp_and_c_wrapper.c -I $(CC_LOC)/../src/ 	

LINKS = chimesFF.o test_cpp_and_c_wrapper.o wrapper-C.o
	
test_cpp_and_c_wrapper : $(LINKS) 
	$(CXX)  $(LINKS) -o test_cpp_and_c_wrapper
	
all:
	make chimesFF.o
	make wrapper-C.o
	make test_cpp_and_c_wrapper.o
	make test_cpp_and_c_wrapper
	
clean:
	rm -f *.o 
	
clean-all:
	rm -f *.o 
	rm -f test_cpp_and_c_wrapper	
