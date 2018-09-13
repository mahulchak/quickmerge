#to compile the executable for quickmerge,the assembly merger

CC = g++
CFLAGS = -g -Wall -std=c++0x -O3 -I .

default: quickmerge
quickmerge: qmergelib.o parser.o errorlib.o quickmerge.o
	$(CC) $(CFLAGS) -o quickmerge qmergelib.o errorlib.o parser.o quickmerge.o

parser.o: parser.cpp qmerge.h
	$(CC) $(CFLAGS) -c parser.cpp

errorlib.o: errorlib.cpp qmerge.h
	$(CC) $(CFLAGS) -c errorlib.cpp

qmergelib.o: qmergelib.cpp qmerge.h
	$(CC) $(CFLAGS) -c qmergelib.cpp

quickmerge.o: quickmerge.cpp qmerge.h
	$(CC) $(CFLAGS) -c quickmerge.cpp

	ln -s merger/quickmerge ../quickmerge
clean:
	$(RM) *.o 
	rm -f quickmerge
	rm -f ../quickmerge
