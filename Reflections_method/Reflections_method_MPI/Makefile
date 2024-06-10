CC = mpic++
CFLAGS = -fsanitize=address -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wno-suggest-attribute=format -lm -ldl -O3

all: main

main: 1.o tr.o cr.o nev.o
	$(CC) $(CFLAGS) 1.o tr.o cr.o nev.o -o main

1.o: 1.cpp
	$(CC) $(CFLAGS) -c 1.cpp

tr.o: tr.cpp
	$(CC) $(CFLAGS) -c tr.cpp

cr.o: cr.cpp
	$(CC) $(CFLAGS) -c cr.cpp

nev.o: nev.cpp
	$(CC) $(CFLAGS) -c nev.cpp

clean:
	rm -rf *.o main
