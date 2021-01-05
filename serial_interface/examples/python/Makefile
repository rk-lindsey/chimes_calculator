C_LOC = $(realpath .)

CXX     = g++ -O3 -std=c++11  # -fPIC

CHIMESFF_LOC=$(C_LOC)/../../../chimesFF/src
CHIMESFF_SRC=$(CHIMESFF_LOC)/chimesFF.cpp
CHIMESFF_HDR=$(CHIMESFF_LOC)/chimesFF.h

chimesFF.o : $(CHIMESFF_SRC)
	$(CXX) -c $(CHIMESFF_SRC) -I $(CHIMESFF_LOC)
	
SERIAL_LOC=$(C_LOC)/../../src
SERIAL_SRC=$(SERIAL_LOC)/serial_chimes_interface.cpp
SERIAL_HDR=$(SERIAL_LOC)/serial_chimes_interface.h

serial_chimes_interface.o : $(SERIAL_SRC) 
	$(CXX) -c $(SERIAL_SRC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC) 

WRAPPER_LOC=$(C_LOC)/../../api
WRAPPER_SRC=$(WRAPPER_LOC)/wrapper-C.cpp
WRAPPER_HDR=$(WRAPPER_LOC)/wrapper-h.cpp

wrapper-C.o : $(WRAPPER_SRC) 
	$(CXX) -c $(WRAPPER_SRC) -I $(WRAPPER_LOC) -I $(SERIAL_LOC) -I $(CHIMESFF_LOC)	
	
	
wrapper-C.so: 	
	$(CXX) -shared -o lib-C_wrapper-serial_interface.so wrapper-C.o serial_chimes_interface.o chimesFF.o	
	
clean:
	rm -f *.o

clean-all:

	make clean
	rm -f *.so

all:
	make chimesFF.o
	make serial_chimes_interface.o
	make wrapper-C.o
	make wrapper-C.so
	make clean
