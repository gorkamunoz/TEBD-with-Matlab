function D = contraction(A, B, C)

D = ncon({A, B, C},{[1 2], [3 1 -1], [3 2 -2]}, [1 2 3]);
