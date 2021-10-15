all:
	make -C build -j`nproc`

test: all
	make -C build test

tags: *.cpp *.h
	ctags --c++-kinds=+p --fields=+iaS --extra=+q --language-force=C++ *

.PHONY: clean tags test

clean:
	make -C build clean
