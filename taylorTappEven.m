function[amp] = taylorTappEven(N,SL,n)
%for this program refer IEEE transation on Antennas and Propagation
%vol. AP-32, no. 10, october 1984. pg 1089-1093
%Taylor patterns for discrete arrays
%by-Alfred T. Villeneuve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta = 10^(SL/20);
Ele = 2*N;
mu0 = cosh((1/(Ele-1))*log((eta + sqrt(eta*eta-1))));
sigma = (n*pi)/((Ele)*acos((1/mu0)*cos((2*n-1)*(pi/(2*(Ele-1))))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cnt = 1:(Ele-1);
    sai_p(cnt) = 2*acos((1/mu0)*cos((2*cnt-1)*(pi/(2*(Ele-1)))));
end
d_sai_p = sigma*sai_p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%numE0
mul = 1;
for q = 1:(n-1)
    mul = mul*(sin(d_sai_p(q)/2))*(sin(d_sai_p(q)/2));
end
numE0 = (Ele)*mul;

%%%%%%%denE0
mul = 1;
for q = 1:(n-1)
    mul = mul*sin((q*pi)/((Ele)))*sin((q*pi)/((Ele)));
end
denE0 = mul;
E0 = numE0/denE0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%E(m2pi/2N)
for m = 1:(n-1)
    num_mul = 1;
    for q = 1:(n-1)
        num_mul = num_mul*sin(((m*pi)/((Ele)))-(d_sai_p(q)/2))*sin(((m*pi)/((Ele)))+(d_sai_p(q)/2));
    end
    
    den_mul =1;
    for q = 1:(n-1)
        if q == m;
            den_mul = den_mul;
        else
            den_mul = den_mul*sin(((m-q)*pi)/((Ele)))*sin(((m+q)*pi)/((Ele)));
        end
    end
    E(m) = (((Ele))*(-1)^m*num_mul)/(sin((m*pi)/((Ele)))*sin((2*m*pi)/((Ele)))*den_mul);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:N
    sum = 0;
    for m = 1:(n-1)
%         sum = sum + E(m)*cos((p*m*2*pi)/((Ele)+1));
        sum = sum + E(m)*cos(((p-0.5)*m*pi)/((1*N)));
    end
%     Amp(p) = (1/((Ele)+1))*(E0 + 2*sum);
    Amp(p) = (1/((Ele)))*(E0 + 2*sum);
end
amp = [fliplr(Amp) Amp];
amp = amp/max(amp);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;
    
 
