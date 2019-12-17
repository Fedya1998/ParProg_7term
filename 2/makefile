x=100
y=100
all: clean old new run

diff: clean new old run dif

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
	mpirun -n 8 2a $(x) $(y)
run_old:
	./2a_old $(x) $(y)

dif:
	diff result_old.txt result.txt
