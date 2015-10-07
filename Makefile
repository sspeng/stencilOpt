size=100
input=jacobi7
opt=_opt

CC=gcc -O3 -fopenmp

time:
	$(CC) $(input)$(opt).c $(input)_timer.c -D SIZE=$(size) 
	./a.out

test:
	$(CC) $(input)$(opt).c $(input)_tester.c -D SIZE=$(size) 
	./a.out

clean:
	rm log_0 log_1 a.out
