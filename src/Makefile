CC=g++

CFLAGS=-c -Wall

src=$(wildcard *.cpp)

obj=$(src: %.cpp=%.o)

target=FVM

$(target):$(obj)
	$(CC) $(obj) -o $(target)

%.o: %.cpp
	$(CC) $< -o $@	

.PHONY: clean

clean:
	rm -rf $(obj) $(target)
