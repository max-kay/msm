main: main.cpp v3.cpp builder.cpp
	g++ -std=c++17 main.cpp -o main -Wall -Wextra -O3 -ffast-math -flto -march=native
