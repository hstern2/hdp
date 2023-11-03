CFLAGS = -g -Wall
CC = gcc

all: test

test: hdp inputFile.fasta out.csv tmp.csv
	diff tmp.csv out.csv

tmp.csv: hdp
	./hdp inputFile.fasta 1-54,272-566 > tmp.csv

hdp: hdp.o util.o hdp_main.o

clean:
	rm -f *.o hdp tmp.csv