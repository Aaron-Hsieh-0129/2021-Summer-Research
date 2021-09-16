# CC and CFLAGS are varilables
CC = g++
CFLAGS = -c
# -c option ask g++ to compile the source files, but do not link.
OPTFLAGS = -O3

all	: bin/vvm2d
	@echo -n ""

bin/vvm2d	: main.o Declare.o Init.o Iteration.o Outputfile.o
			$(CC) $(OPTFLAGS) main.o Declare.o Init.o Outputfile.o Iteration.o -o bin/vvm2d
main.o 	   	: src/main.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Declare.o	: src/Declare.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Init.o		: src/Init.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Outputfile.o: src/Outputfile.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
Iteration.o	: src/Iteration.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@


# clean all the .o and executable files
clean:
		rm -rf *.o bin/*
		rm -rf outputs/init/* outputs/qc/* outputs/qr/* outputs/qv/* outputs/th/* outputs/th2d/* outputs/u/* outputs/w/* outputs/zeta/*
		rm -rf graphs/init/* graphs/qc/* outputs/qc+qr/* outputs/qc+qr+th+u+w/* outputs/qr+th+u+w/* outputs/qv+qc/* graphs/qr/* graphs/qv/* graphs/th/* graphs/th2d/* graphs/u/* graphs/w/* graphs/zeta/*