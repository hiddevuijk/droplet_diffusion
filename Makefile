
TARGET = test.exe
OBJS = main.o 
CC = g++

# boost include for mac
INC=-I/opt/homebrew/Cellar/boost/1.80.0/include

#CFLAGS = -c -Wall  -g -std=c++11
#LFLAGS = -Wall  -g -std=c++11
CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
LFALGS = -Wall -O3 -DNDEBUG -std=c++11

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

main.o: main.cpp droplet_diffusion.h concentration.h config_file.h
	$(CC) $(CFLAGS) $(INC)  main.cpp


.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: cleandata
cleandata:
	rm -f *.dat data/*
  
