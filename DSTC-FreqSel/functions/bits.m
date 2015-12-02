% random bits generation

function out=bits(N)
    out=rand(N,1)>0.5;
end