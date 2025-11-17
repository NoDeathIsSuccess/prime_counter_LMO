# prime_counter_LMO
A LMO-like algorithm.

The file `primecounter.cpp` contains two parts: calculating $π(x)$ and calculating the $10$ largest prime numbers less than or equal to $x$.

The method used to calculate $π(x)$ is a kind of combinatorial prime counting algorithms, which is programmed in `lmo_pi.h`.

The $10$ largest prime numbers less than or equal to $x$ are obtained by simply traversing backwards and performing Miller-Rabin primality tests on each one.

`algcmp1.h` includes another algorithm implementation, but it is slightly less efficient.

`verifier.h` uses a simple sieve method, which is only used to verify the correctness of the algorithm in small-scale cases $(\le 10^{6})$.
