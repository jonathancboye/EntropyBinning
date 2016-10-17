all: main.o gene.o
	g++ main.o gene.o -o entropyCRsfs.exe 
main.o: main.cpp gene.h
	g++ -c main.cpp 
gene.o: gene.cpp gene.h
	g++ -c gene.cpp
clean:
	rm -f *~ *# *.o
