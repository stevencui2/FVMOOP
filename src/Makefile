CC=g++

CFLAGS=-c -Wall

output: main.o forAllOperations.o Grid.o Fields.o Solution.o
	$(CC) main.o forAllOperations.o Grid.o Fields.o Solution.o -o FVM
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
forAllOperatons.o: forAllOperations.cpp forAllOperations.h
	$(CC) $(CFLAGS) forAllOperations.cpp forAllOperations.h
Grid.o: Grid.cpp Grid.h
	$(CC) $(CFLAGS) Grid.cpp Grid.h
Fields.o: Fields.cpp Fields.h
	$(CC) $(CFLAGS) Fields.cpp Fields.h
Solution.o: Solution.cpp Solution.h
	$(CC) $(CFLAGS) Solution.cpp Solution.h
clean:
	rm -rf *.o *.gch FVM
