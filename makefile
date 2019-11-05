CXX = g++
CC  = $(CXX)

CXXFLAGS = -std=c++1z

# Partie commentÃ©e : choisissez les options que vous voulez avoir
#                    en dÃ©commentant la/les lignes correspondantes
#
CXXFLAGS += -pedantic -Wall #       # pour les purs et durs
CXXFLAGS += -g #                    # pour debugger
# CXXFLAGS += -pg #                 # pour profiler
# LDFLAGS  += -pg #                 # pour profiler
# CXXFLAGS += -O2 #                 # pour optimiser la vitesse

LIB=LIB/
SRC=SRC/

all: main

$(LIB)main.o: $(SRC)main.cc
	$(CC) $(CXXFLAGS) -c -o $@ $<

#$(LIB)Vecteur.o: $(SRC)Vecteur.cc $(SRC)Vecteur.h $(SRC)Erreur.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<
#
#$(LIB)Systeme.o : $(SRC)Systeme.cc $(SRC)Systeme.h $(SRC)Vecteur.h $(SRC)Erreur.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<
#
#$(LIB)Exercice4.o : $(SRC)Exercice4.cc $(SRC)Systeme.h $(SRC)Vecteur.h $(SRC)Erreur.h
#	$(CC) $(CXXFLAGS) -c -o $@ $<
#
#$(LIB)performance.o : $(SRC)performance.cc
#	$(CC) $(CXXFLAGS) -c -o $@ $<

#Exercice4: $(LIB)Exercice4.o $(LIB)Systeme.o $(LIB)Vecteur.o $(LIB)Erreur.o
#	$(CC) $(CXXFLAGS) $(LIB)Exercice4.o $(LIB)Systeme.o $(LIB)Vecteur.o $(LIB)Erreur.o -o Exercice4

main: $(LIB)main.o
	$(CC) $(CXXFLAGS) $(LIB)main.o -o main
