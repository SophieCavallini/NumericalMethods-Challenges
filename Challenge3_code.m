clear all
close all
clc
format long
%% part 1
R = 3;
V = 10;
A = [3*R -R -R 0;
    -R +3*R 0 -R;
    -R 0 3*R 0;
    0 -R 0 3*R];
b = [V -V V -V]';
% System to be solved: A*x=b
x = A\b;
x = inv(A)*b;
x = A^(-1)*b;
[L, U, P] = lu(A);
y = fwsub(L, b);
x = bksub(U, y);
%% part 2
clear all
R = 3;
r = @(i) R^(-i);
V = 10;
v = @(i) V^(-i+1);
A = [3*r(0) -r(0) 0 0 -r(0) 0 0 0;
    -r(0) 3*r(0) 0 0 0 -r(0) 0 0;
    0 0 3*r(1) -r(1) 0 0 -r(1) 0;
    0 0 -r(1) 3*r(1) 0  0 0 -r(1);
    -r(0) 0 0 0 3*r(0) 0 0 0;
    0 -r(0) 0 0 0 2*r(0)+r(1) -r(1) 0;
    0 0 -r(1) 0 0 -r(1) 2*r(1) 0;
    0 0 0 -r(1) 0 0 0 2*r(1)+r(2)];
b = [v(0) -v(1) v(1) -v(1) v(0) -v(0) v(1) -v(1)]';
% System to be solved Ax = b
x = A\b
x = inv(A)*b
x = A^(-1)*b
[L, U, P] = lu(A);
y = fwsub(L, b);
x = bksub(U, y)
%% part 3
[A3, b3] = buildAb(3);
[A5, b5] = buildAb(5);
[A7, b7] = buildAb(7);
x3 = A3\b3
x5 = A5\b5
x7 = A7\b7

%% functions

function x = bksub(U, b)
if (U ~= triu(U))
    disp('Error: L is not upper triangular');
    x = NaN;
    return;
end
n = length(b);
x = zeros(n, 1);
x(n) = b(n)/U(n,n);
for i = 1:(n-1)
    k = n - i;
    somma = 0;
    for s = k+1:n
        somma = somma + U(k,s)*x(s);
    end
    x(k)=(b(k)-somma)/U(k,k);
end
end

function x = fwsub(L,b)
if (L ~= tril(L))
    disp('Error: L is not lower triangular');
    x = NaN;
    return;
end
n = length(b);
x = zeros(n, 1);
x(1) = b(1)/L(1,1);
for k = 2:n
    somma = 0;
    for s = 1:k-1
        somma = somma + L(k,s)*x(s);
    end
    x(k)=(b(k)-somma)/L(k,k);
end
end

function aa = block4A(i)
R = 3;
r = @(j) R^(-j);
aa = zeros(4,6);
aa = [0 3*r(i) -r(i) -r(i) 0 0;
      0 -r(i) 3*r(i) 0 -r(i) 0;
      -r(i) -r(i) 0 3*r(i) 0 0;
      0 0 -r(i) 0 2*r(i+1)+r(i) r(i+1)];
end

function [A, b] = buildAb(N)
V = 10;
v = @(j) V^(-j+1);
b = zeros(4*(N+1),1);
A = zeros(4*(N+1),4*(N+1)+2);
for i=1:2*(N+1)
    if mod(i,4) == 2
        b(i) = -v(floor((i-1)/4)+1);
    else
        b(i) = v(floor((i-1)/4))*(-1)^(mod(i-1,2));
    end
end
for i=0:N
    A(4*i+1:4*(i+1),4*i+1:4*(i+1)+2)=block4A(i);
end
A = A(:,2:4*(N+1)+1);
end
