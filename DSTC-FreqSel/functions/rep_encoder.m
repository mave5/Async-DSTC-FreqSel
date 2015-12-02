% repetition encoder

function out=rep_encoder(b_in,r)
% b_in : binary input data
% r    : code rate
% out  : coded data

out=reshape(repmat(b_in,r,1),[],1);

end



