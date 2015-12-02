function alfa=raised_cosine(tau,a)
% tau: delay
% a : roll-off factor
    alfa=sinc(tau).*cos(pi*a*tau)./(1-4*a^2*tau^2);
end