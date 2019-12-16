CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o blocchi.o GA.o
	$(CC) random.o blocchi.o GA.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
GA.o: GA.cpp GA.h random.cpp random.h
	$(CC) -c GA.cpp -o GA.o $(CFLAGS)
blocchi.o : blocchi.cpp blocchi.h
	$(CC) -c blocchi.cpp -o blocchi.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
