makemandelbrot: main.cpp
	g++ main.cpp -O3 -DNDEBUG -o mandelbrot -lncurses
