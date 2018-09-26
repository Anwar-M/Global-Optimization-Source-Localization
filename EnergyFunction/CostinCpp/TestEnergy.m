function z = TestEnergy(Xs,cpreal,cpimag,f,csound,XA,YA,ZA)

% The diagonal are the powers! So amplitude squared. To obtain amplitude
% take square root in script.

% Two sources
X = [Xs(1:4); Xs(5:8)];
% Four sources
% X = [Xs(1:4); Xs(5:8); Xs(9:12); Xs(13:16)];

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
format long;
disp(C_mod);
format short;
% z = sum(sum( (real(C_mod) - squeeze(cpreal(fi,:,:))).^2 )) + ...
%     sum(sum( (imag(C_mod) - squeeze(cpimag(fi,:,:))).^2 ));

z = sum(sum( (real(C_mod) - cpreal).^2 )) + ...
    sum(sum( (imag(C_mod) - cpimag).^2 ));
