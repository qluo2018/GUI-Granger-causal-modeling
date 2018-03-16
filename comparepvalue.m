%%
function pvalue = comparepvalue(s, p, Np, q, Nq)
sigma = p*(1-p)/Np + q*(1-q)/Nq;
pvalue = normcdf(-abs(s)/sqrt(sigma));