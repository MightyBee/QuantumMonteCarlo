CXX = g++
CC  = $(CXX)

CXXFLAGS = -std=c++1z

# Partie commentÃ©e : choisissez les options que vous voulez avoir
#                    en dÃ©commentant la/les lignes correspondantes
#
CXXFLAGS += -pedantic -Wall -Wextra#       # pour les purs et durs
CXXFLAGS += -g #                    # pour debugger
# CXXFLAGS += -pg #                 # pour profiler
# LDFLAGS  += -pg #                 # pour profiler
# CXXFLAGS += -O2 #                 # pour optimiser la vitesse

LIB=LIB/
SRC=SRC/
BIN=BIN/

all: $(BIN)main2 #metropolis #main

#$(LIB)main.o : $(SRC)main.cc $(SRC)ConfigFile.tcc $(SRC)ConfigFile.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<

$(LIB)main2.o : $(SRC)main2.cc
	$(CC) $(CXXFLAGS) -c -o $@ $<

#$(LIB)ConfigFile.o : $(SRC)ConfigFile.tcc $(SRC)ConfigFile.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<

#$(LIB)potential.o : $(SRC)potential.cc $(SRC)potential.h $(SRC)instance.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<

#$(LIB)system.o : $(SRC)system.cc $(SRC)system.h $(SRC)potential.h $(SRC)instance.h #$(SRC)ConfigFile.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<

#$(LIB)metropolis.o : $(SRC)metropolis.cc $(SRC)potential.h $(SRC)system.h $(SRC)instance.h #$(SRC)ConfigFile.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<

#metropolis: $(LIB)metropolis.o #$(LIB)system.o $(LIB)potential.o
#	$(CC) $(CXXFLAGS) $(LIB)metropolis.o $(LIB)system.o $(LIB)potential.o -o metropolis

#main: $(LIB)main.o
#	$(CC) $(CXXFLAGS) $(LIB)main.o -o main

$(BIN)main2: $(LIB)main2.o
	$(CC) $(CXXFLAGS) $(LIB)main2.o -o $(BIN)main2
