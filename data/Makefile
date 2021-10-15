NAME?=convert
CXX=clang++
LDFLAGS+= -lpthread -lomp -lfftw3f -lfftw3 -lHadrons -lGrid -lcrypto -llime -lhdf5_cpp -lz -lmpi -lhdf5
CPPFLAGS+= -g -Wall -std=c++11 -O2 -march=native -I./tclap/include

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	LDFLAGS+= -lrt -lmpi_cxx -lhdf5
	CPPFLAGS+= -fopenmp
else
	LDFLAGS+=
endif

CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)

all: $(NAME)

$(NAME): $(OBJ_FILES)
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp *.hpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<

tags: *.cpp *.h
	ctags --c++-kinds=+p --fields=+iaS --extra=+q --language-force=C++ -R .

.PHONY: clean tags

clean:
	rm -f *.o $(NAME)
