CC=g++ 
CFLAGS= -w -Wall -I ./eigen/ -std=c++11
all:
	$(CC) $(CFLAGS) localization_v2.4_dynamic.cpp  -o locrun

.PHONY: opt
opt:override CFLAGS += -O2
opt:all

.PHONY: prof
prof:override CFLAGS += -p
prof:all

.PHONY: clean
clean:
	rm locrun
