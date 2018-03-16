function pvalue = sumF2(x, n, d1, d22)
% significance of the summation of n F distributions
for t = 1 : 1000
    sample(t) = 0;
    for i = 1 : n
        sample(t) = sample(t) + frnd(d1,d22(i));
    end
end
PD = fitdist(sample', 'kernel');
pvalue = 1-cdf(PD, x);