all: task1 task2
	rm -rf Matrix/matrix.o

task2: task2.o matrix.o
	g++ task2.o Matrix/matrix.o -o task2
	rm -rf task2.o

task1: task1.cpp
	g++ -std=c++17 -Wall -pedantic task1.cpp -o task1

task2.o: task2.cpp
	g++ -std=c++17 -Wall -pedantic -c task2.cpp

matrix.o: Matrix/matrix.cpp
	g++ -std=c++17 -Wall -pedantic -c Matrix/matrix.cpp -o Matrix/matrix.o


clean:
	rm -rf *.o task1 task2 Matrix/matrix.o .ipynb_checkpoints