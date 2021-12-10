% To get the answer for Problem 4(a) and 4(b)
clear all;close all;clc;
%(a)
beta_CDF(0.25,3,3)
beta_CDF(0.5,4,5)

%(b)
alpha = 5;
beta = 6;
mean = alpha/(alpha+beta);
j = 1;
error = zeros(5,1);
for N = [10 50 100 500 1000]
x = Sampling_beta(N,alpha,beta);
sample_mean = sum(x)/N;
error(j) = abs(sample_mean-mean);
j = j+1;
clear x;
end
%
N = [10 50 100 500 1000];
plot(N,error,'ko--');
xlabel('N');
ylabel('absolute error between sample mean and mean');
