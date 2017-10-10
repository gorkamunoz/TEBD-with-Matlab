function mps = initial_mps(N,d,D)

mps = cell(1,N);
for i = 1:N/2
    mps{i} = rand(d,min([d^(i-1), D]),min([d^i, D]))+1i*rand(d,min([d^(i-1), D]),min([d^i, D]));
    mps{N+1-i} = rand(d,min([d^i, D]),min([d^(i-1), D]))+1i*rand(d,min([d^i, D]),min([d^(i-1), D]));
end


%% Normalizing the MPS

A = mps{N};
[q, r] = qr(A,0);
mps{N} = q;
mps{N-1} = ncon({mps{N-1},ctranspose(r)}, {[-1 -2 1], [1 -3]});

for i=N-1:-1:2
    
    A = mps{i};
    shape = size(permute(A, [2 1 3]));
    B = reshape(permute(A, [2 1 3]), shape(1), []);
    [q, r] = qr(ctranspose(B),0);
    mps{i} = permute(reshape(ctranspose(q),shape), [2 1 3]);
    mps{i-1} = ncon({mps{i-1},ctranspose(r)}, {[-1 -2 1], [1 -3]});
    
end

A = mps{1};
shape = size(A);
B = reshape(A, 1, []);
[q, ~] = qr(ctranspose(B),0);
mps{1} = reshape(ctranspose(q),shape);

%ncon({mps{5},permute(mps{5},[1 3 2])},{[1 -1 2], [1 2 -2]})


end