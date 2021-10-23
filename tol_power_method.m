function [lambda] = tol_power_method(A)
%Initialize
w = randn(size(A,1),1);
tol = 1e-4;
v_new = ones(size(A,1),1);
while 1
    w = A*w;
    v = w / norm(w);
    if norm(v-v_new) > tol
        break
    end
       v_new = v;

end
    lambda = (v'*A*v) / (v'*v);
    fprintf("maximum eigenvalue %4f", lambda)