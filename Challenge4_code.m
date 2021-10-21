clear all
close all
clc
format long

% CHALLENGE 3: generate multivariate Gaussian data

% %Generate N realizations of an independent d􀀀dimensional unit normal
% random vector using MATLAB. Transform them into N realizations of a
% d-dimensional correlated Gaussian random vector with prescribed mean
% vector mu and a covariance matrix Sigma: 
% Compute the empirical mean vector mu_emp and covariance matrix Sigma_ emp
% and compute the Mahalanobis distance of a generic d-dimensional random vector 
% from the empirical mean vector.
N = 500; %number of realizations
d = 2; %dimension of the random vector


A = randn(N, d); %independent unit normal random vector of dimension d 
                 %with N realizations
mu = 9.4+5.4*randn(1,d); %prescribed mean vector
bigsigma = [1 .5; .5 2]; %prescibed covariance vector
L = chol(bigsigma);
Z = repmat(mu, N, 1) + randn(N, d)*L; %correlated Gaussian random vector
mu_empZ = sum(Z)/N; %empirical mean of data
Zcent = zeros(N,d);
%center the global data matrix around the mean
for i=1:N
    Zcent(i,:)=Z(i,:)-mu_empZ;
end
bigsigma_empZ=Zcent'*Zcent/(N-1);%empirical covariance matrix
L_empZ = chol(bigsigma_empZ,'lower');
k = 3; %data sample that will be perturbed
z = Z(k,:)' + 4*randn(size(d,1));
Mahalanobis_distance=sqrt((z-mu_empZ')'*(L_empZ\(z-mu_empZ'))); %compute the Mahalanobis distance
fprintf('Mahalanobis distance: %f ', Mahalanobis_distance);