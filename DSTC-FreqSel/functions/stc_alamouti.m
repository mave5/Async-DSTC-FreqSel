% space time code, Alamouti code
% [u1 -u2';u2 u1']
% in: input symbols, any constelation such as QPSK symbols
% out: Alamouti codewords,matrics of 2*2

function out=stc_alamouti(in,type)

if nargin==1
    type=1;
end

u1=in(1:2:end);%odd symbols
u2=in(2:2:end);%even symbols

u{length(u1)}=[];% define an empty cell, memory of codewords

% creat the codewords
for k=1:length(u1)
    if type==1
    u{k}=[u1(k) -conj(u2(k));u2(k) conj(u1(k))]/(sqrt(abs(u1(k))^2+abs(u2(k))^2));
    else
    u{k}=[u1(k) (u2(k));-conj(u2(k)) conj(u1(k))]/(sqrt(abs(u1(k))^2+abs(u2(k))^2));
    end
end
out=u;
end


% type=1
% [s1 -s2']
% [s2  s1']

% type=2
% [s1    s2]
% [-s2'  s1']