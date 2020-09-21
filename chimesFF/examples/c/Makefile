CC_LOC = $(realpath .)

CXX=g++ -O3 -std=c++11 

CHIMESFF_SRC=$(CC_LOC)/../../src/chimesFF.cpp
CHIMESFF_HDR=$(CC_LOC)/../../src/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC)
		
WRAPPER_SRC=$(CC_LOC)/../../api/wrapper-C.cpp
WRAPPER_HDR=$(CC_LOC)/../../api/wrapper-C.h
		
wrapper-C.o : $(WRAPPER_SRC) $(WRAPPER_HDR) $(CHIMESFF_SRC)  $(CHIMESFF_HDR)
	$(CXX) -c $(WRAPPER_SRC) $(CHIMESFF_SRC) -I $(CC_LOC)/../../api/ -I $(CC_LOC)/../../src/
	
CC = gcc

test_wrapper-C.o : main.c $(WRAPPER_SRC) $(WRAPPER_HDR)
	$(CC) -c main.c -o test_wrapper-C.o -I $(CC_LOC)/../../api/ 	

LINKS = chimesFF.o test_wrapper-C.o wrapper-C.o
	
test_wrapper-C : $(LINKS) 
	$(CXX)  $(LINKS) -o test_wrapper-C
	
clean:
	rm -f *.o 
	
clean-all:
	rm -f *.o 
	rm -f test_wrapper-C

all:
	make chimesFF.o
	make wrapper-C.o
	make test_wrapper-C.o
	make test_wrapper-C
	make clean
