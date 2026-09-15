all:
	mkdir -p ./build
	gcc -O3 -march=native -I/usr/local/include/suitesparse ./src/masked_states.c -o ./build/masked -lgraphblas -llagraph -llagraphx
	gcc -O3 -march=native -I/usr/local/include/suitesparse ./src/straightforward.c -o ./build/straight -lgraphblas -llagraph -llagraphx
