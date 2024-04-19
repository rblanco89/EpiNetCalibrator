COMPILER = g++
LIBS = -lm
CC = $(COMPILER) 
OBJS = mainEpiNet.o

EpiNetSimulator : $(OBJS)
	$(CC) $^ -o $@ $(LIBS)

%.o : %.cpp
	$(CC) -c $< -o $@

.PHONY: clean

clean:
	rm -f *.o *~
