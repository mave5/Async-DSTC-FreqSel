% differenrial encoder MIMO-OFDM

function s_k=diff_encoder_v(V_in,s_km1)

N=length(V_in);

for t=1:N
    s_k(:,t)=V_in{t}*s_km1(:,t);
end

end


