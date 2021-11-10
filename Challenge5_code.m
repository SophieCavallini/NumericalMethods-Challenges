clear all
close all
clc
format long
clear all

[A25, b25] = buildAb(25);
x25 = A25\b25;

apriori = cond(A25)*eps;
res25 = b25-A25*x25;
apost = cond(A25)*norm(res25,2)/norm(b25,2);

i = 25;

while apost < 0.1
    i = i + 1;
    [A, b] = buildAb(i);
    x = A\b;
    res = b - A*x;
    apost = cond(A)*norm(res,2)/norm(b,2);
end
fprintf('Maximum number of loop blocks: %d\n',i);

A = buildAbp(25);
apriori = cond(A)*eps


%% functions
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
v = @(j) V^((-j+1)/100)  + 0.01*V^((-j+1)/100)*rand(1);
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

function aa = block4Ap(i)
R = 3;
r = @(j) R^(-j)*(1+0.01*rand(1));
aa = zeros(4,6);
aa = [0 3*r(i) -r(i) -r(i) 0 0;
      0 -r(i) 3*r(i) 0 -r(i) 0;
      -r(i) -r(i) 0 3*r(i) 0 0;
      0 0 -r(i) 0 2*r(i+1)+r(i) r(i+1)];
end

function [A, b] = buildAbp(N)
V = 10;
v = @(j) V^((-j+1)/100)  + 0.01*V^((-j+1)/100)*rand(1);
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
    A(4*i+1:4*(i+1),4*i+1:4*(i+1)+2)=block4Ap(i);
end
A = A(:,2:4*(N+1)+1);
end
