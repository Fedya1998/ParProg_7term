all: clean new old run diff

nodiff: clean new old run

old: 2a_old.c
	mpicc 2a_old.c -lm -std=c99 -o 2a_old

new: 2a.c
	mpicc 2a.c -lm -std=c99 -o 2a

clean:
	rm -rf 2a 2a_old *.txt
clean_old:
	rm -rf 2a_old result_old.txt
clean_new:
	rm -rf 2a result.txt

run: run_new run_old

run_new:
	mpirun -n 4 2a 100 100
run_old:
	./2a_old

diff:
	diff result_old.txt result.txt
