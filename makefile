CC := gcc
override CFLAGS += --std=gnu89 -Wall -Wextra -pedantic -O3
LDFLAGS := -lm

EXEC := pagerank hits

all: $(EXEC)

compdb:
	bear -- make clean all

pagerank: pagerank.o utils.o
	$(CC) -o pagerank pagerank.o utils.o $(CFLAGS) $(LDFLAGS)

hits: hits.o jaccard.o utils.o
	$(CC) -o hits hits.o jaccard.o utils.o $(CFLAGS) $(LDFLAGS)

pagerank.o: src/pagerank.c src/utils.h
	$(CC) -c src/pagerank.c $(CFLAGS)

hits.o: src/hits.c src/jaccard.h src/utils.h
	$(CC) -c src/hits.c $(CFLAGS)

utils.o: src/utils.c src/utils.h
	$(CC) -c src/utils.c $(CFLAGS)

jaccard.o: src/jaccard.c src/jaccard.h
	$(CC) -c src/jaccard.c $(CFLAGS)

clean:
	rm -rf *.o $(EXEC) *.pr *.hits HITS_* PR_* *.csv
