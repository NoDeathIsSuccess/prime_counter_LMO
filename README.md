# prime_counter_LMO
A lightly LMO-like algorithm.

The file `primecounter.cpp` contains two parts: calculating $π(x)$ and calculating the $10$ largest prime numbers less than or equal to $x$.

The method used to calculate $π(x)$ is a kind of combinatorial prime counting algorithms, which is programmed in `lmo_pi.h`.

The 10 largest prime numbers less than or equal to $x$ are obtained by simply traversing backwards and performing Miller-Rabin primality tests on each one.

`algcmp1.h` includes another algorithm implementation, but it is slightly less efficient.

`verifier.h` uses a simple sieve method, which is only used to verify the correctness of the algorithm in small-scale cases ($\le 10^{6}$).

No additional libraries. No detailed optimizations have been performed. OpenMP parallelism is used simply in one place, and the improvement is very limited. Perhaps performance could be further improved.

# Benchmark
| $x$         | $\pi(x)$      | Time (ms)   |
| ----------- | ------------- | ----------- |
| 1e8         | 5761455       | 4.03        |
| 1e9         | 50847534      | 5.40        |
| 1e10        | 455052511     | 7.20        |
| 1e11        | 4118054813    | 17.09       |
| 1e12        | 37607912018   | 55.68       |
| 1e13        | 346065536839  | 255.04      |
| 1e14        | 3204941750802 | 1334.72     |

Intel Core i5-12600KF Processor, up to 4.90 GHz.


