
#THIS IS A MAKEFILE FOR HEROKU ONLY, local_make has the one for this machine

CC = g++ -I /app/vendor/vips/include
CCOPTS = -g
LINK = g++ -L /app/vendor/vips/lib

rscft: rscft.o main.o
	$(LINK) SCFT_real/rscft.o SCFT_real/main.o -lfftw3 -lm -o rscft
clean:
	-rm SCFT_real/rscft.o SCFT_real/main.o SCFT_real/rscft
rscft.o: rscft.cc rscft.h
	$(CC) $(CCOPTS) -c SCFT_real/rscft.cc
main.o: main.cc rscft.h
	$(CC) $(CCOPTS) -c SCFT_real/main.cc
