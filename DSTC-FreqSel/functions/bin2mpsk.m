% this function modulates binary data to M-psk symbols
% gray mapping 

% M=2 BPSk
% 1 -> 1
% 0 -> -1
% M=4 QPSK
% 00 -> 1
% 01 -> j
% 11 -> -1
% 10 -> -j

% M=8, 8-PSK
% 000 -> -1-j
% 001 -> -1
% 010 -> j
% 011 -> -1+j
% 100 -> -j
% 101 -> 1-j
% 110 -> 1+j
% 111 -> 1


function s= bin2mpsk(bits,M,rot_ang) 
% bits: binary 0/1 input data
% M   : MPSK
% s   : MPSK symbols exp(j*2*pi*m/M)
% rot_ang: rotation angle
if nargin==2
    rot_ang=0;
end

    if M==2
        s(:,1)=(2*bits-1)*exp(1i*rot_ang);
    elseif M==4    
        odd=bits(1:2:end);%odd bits
        even=bits(2:2:end);%even bits
        s=(1-odd-even)+ 1i*(even-odd);
        s=s*exp(1i*rot_ang);
    elseif M==8
        s=zeros(length(bits)/log2(M),1);
        for k=1:3:length(bits)
            l=floor(k/3)+1;
            b2d=bits(k)*4+bits(k+1)*2+bits(k+2);
            switch b2d
                case 0
                    s(l)=(-1-1i)/sqrt(2);
                case 1
                    s(l)=-1;
                case 2
                    s(l)=1i;
                case 3
                    s(l)=(-1+1i)/sqrt(2);
                case 4
                    s(l)=-1i;
                case 5
                    s(l)=(1-1i)/sqrt(2);
                case 6
                    s(l)=(1+1i)/sqrt(2);
                case 7
                    s(l)=1;
            end    
        end
    end
    
end
