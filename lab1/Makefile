all: task1 task2 task3 task4 task5

task1: task1.o matrix.o
	g++ task1.o matrix.o -o task1
	rm -rf task1.o

task2: task2.o matrix.o
	g++ task2.o matrix.o -o task2
	rm -rf task2.o

task3: task3.o matrix.o
	g++ task3.o matrix.o -o task3
	rm -rf task3.o

task4: task4.o matrix.o
	g++ task4.o matrix.o -o task4
	rm -rf task4.o

task5: task5.o matrix.o
	g++ task5.o matrix.o -o task5
	rm -rf task5.o

task1.o: task1.cpp
	g++ -std=c++17 -Wall -pedantic -c task1.cpp

task2.o: task2.cpp
	g++ -std=c++17 -Wall -pedantic -c task2.cpp

task3.o: task3.cpp
	g++ -std=c++17 -Wall -pedantic -c task3.cpp

task4.o: task4.cpp
	g++ -std=c++17 -Wall -pedantic -c task4.cpp

task5.o: task5.cpp
	g++ -std=c++17 -Wall -pedantic -c task5.cpp

matrix.o: matrix.cpp
	g++ -std=c++17 -Wall -pedantic -c matrix.cpp


clean:
	rm -rf *.o task1 task2 task3 task4 task5