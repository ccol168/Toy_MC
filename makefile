CC = g++
CFLAGS = -Wall -O3 --std=c++11

Toy_MC.exe : Toy_MC.cpp
	$(CC) Toy_MC.cpp -o Toy_MC.exe
