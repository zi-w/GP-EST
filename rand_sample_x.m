function initx = rand_sample_x(N, D, xlimit)
% randomly sample N points in the space specified by xlimit
initx = zeros(N,D);
for j = 1:D
    initx(:,j) = rand(N,1)*(xlimit(2,j)-xlimit(1,j)) + xlimit(1,j);
end