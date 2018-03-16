function LLF = arma_LLF(errors, precov, Nr, Nl, order) 
[k,T] = size(errors');
LLF = 0;
likelihoods = zeros(T,1);
for i = 1 : Nr
    for j = Nl(i,1)+order : Nl(i,2)
        likelihoods(j) = k*log(2*pi)+(log(det(precov)) + errors(j,:)*precov^(-1)*errors(j,:)');
        LLF = LLF+likelihoods(j);        
    end                                                                                                                                                                                   
end    
LLF = -0.5*(LLF);