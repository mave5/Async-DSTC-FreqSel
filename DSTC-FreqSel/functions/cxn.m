% Complex Normal generator
% N: numbers
% N0: variance of the Complex Normal
function out=cxn(N,N0)
    out=sqrt(N0/2)*(randn(1,N)+i*randn(1,N));
end