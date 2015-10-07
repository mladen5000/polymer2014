
#THIS IS A MAKEFILE FOR HEROKU ONLY, local_make has the one for this machine

CC = g++ -I /app/vendor/vips/include
CCOPTS = -g
LINK = g++ -L /app/vendor/vips/lib

rscft: rscft.o main.o
	$(LINK) rscft.o main.o -lfftw3 -lm -o rscft
clean:
	-rm rscft.o main.o rscft
rscft.o: rscft.cc rscft.h
	$(CC) $(CCOPTS) -c rscft.cc
main.o: main.cc rscft.h
	$(CC) $(CCOPTS) -c main.cc
