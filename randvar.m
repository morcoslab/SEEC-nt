function X=randvar(P,n)
% returns a realization for a discreate random variable given its
% cumulative probability distribution.
[~,X] = histc(rand(1,n),P);
X=X+1;
end