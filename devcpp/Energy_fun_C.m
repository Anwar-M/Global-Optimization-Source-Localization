function z = Energy_fun_C(Xs,cpreal,cpimag,f,csound,XA,YA,ZA)

% The diagonal are the powers! So amplitude squared. To obtain amplitude
% take square root in script.

% Two sources
X = [Xs(1:4); Xs(5:8)];

n_src = size(X, 1);
n_mic = length(XA);

% Calculate steering matrix
A = zeros(n_mic, n_src);
for i = 1:n_src
    RR = sqrt((XA-X(i,1)).^2 + (YA-X(i,2)).^2 + (ZA-X(i,3)).^2);
    A(:, i) = exp(-2*pi*1i*f.*RR/csound)./(4*pi*RR);
end
D = diag(X(:,4),0);
C_mod = 0.5*A*D*A';

z = sum(sum( (real(C_mod) - cpreal).^2 )) + ...
    sum(sum( (imag(C_mod) - cpimag).^2 ));

% for k1 = 1:length(RR)
%     for k3 = 1:length(RR)
%         C(k1,k3) = 0.5*a^2*exp(-2*pi*i*f*(RR(k3)-RR(k1))/csound);
%         C(k1,k3) = C(k1,k3)/(4*pi*RR(k1)*4*pi*RR(k3));
%     end
% end
% z = sum(sum(abs(imag(C)-squeeze(cpimag(fi,:,:))) + abs(real(C)-squeeze(cpreal(fi,:,:)))));
% z = sum(sum(abs(imag(C)-squeeze(cpimag(fi,:,:))).^2 + abs(real(C)-squeeze(cpreal(fi,:,:))).^2));
