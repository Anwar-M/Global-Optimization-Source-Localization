function z = Energy_fun_B(X,cpreal,cpimag,fi,f,csound,XA,YA,ZA)

XS = X(1);
YS = X(2);
ZS = X(3);
a = X(4);
if (length(X) == 5)
    csound = X(5);
end

% Conform paper of Mexcklenbrauker, JASA 1999: However, will not find the
% source amplitude, see program testje_normalization.m for check on
% normalization
RR = sqrt((XA-XS).^2 + (YA-YS).^2 + (ZA-ZS).^2);
d = (a./(4*pi*RR).*exp(2*pi*1i*f*RR/csound)).';
d = d/norm(d);
z = 1 - d'*(squeeze(cpreal(fi,:,:)) + 1i*squeeze(cpimag(fi,:,:)))*d/abs(trace(squeeze(cpreal(fi,:,:)) + 1i*squeeze(cpimag(fi,:,:))));
if imag(trace(squeeze(cpreal(fi,:,:)) + 1i*squeeze(cpimag(fi,:,:)))) > 0.0001
    keyboard
end
z = real(z);

% End of program