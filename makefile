CXX ?= g++
EXEC= intervalmiss
SRC_DIR := src
OBJ_DIR := obj

CXX ?= g++

CFLAGS = -O3 -std=c++11 -lpthread

SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

.PHONY: test clean deb


all: $(EXEC) 

clean:
	-rm -f  $(EXEC) $(OBJ_FILES)

$(EXEC): $(OBJ_FILES)
	$(CXX)  -o $@ $^ $(CFLAGS)



$(OBJ_DIR):
	mkdir $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(OBJ_DIR)
	$(CXX)   -c -o $@ $< $(CFLAGS)

#to compile the reading test from file
#g++  -O3 -std=c++11  read-file-fast.cc SAMR.cpp  -o read-file-fast
