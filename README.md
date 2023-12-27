# PSO for C/CPP using EIGEN 3.4.0
 
# Algorithm 1: Particle Swarm Optimization 

```plaintext
1: Initialize population of particles with random positions and velocities
2: Initialize global best position vector: P_global
3: while stopping criteria are not met do
4:     for each particle do
5:         Update particle’s best position vector: p_i
6:         Update global best position vector: P_global
7:         for each dimension do
8:             Update particle’s velocity: v_i(t + 1) = w * v_i(t) + c1 * r1 * (p_i(t) - x_i(t)) + c2 * r2 * (P_global(t) - x_i(t))
9:             Update particle’s position: x_i(t + 1) = x_i(t) + v_i(t + 1)
10:         end for
11:     end for
12: end while
13: return global best position: P_global
```
