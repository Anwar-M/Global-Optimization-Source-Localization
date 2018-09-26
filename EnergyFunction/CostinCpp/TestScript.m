clc;
Xs = [1.1, 2.1, 1.0, 1.0, 0.5, 1.1, 1.0, 2.0];
micpos = [0.5, 0.5, 0; 1.0, -0.4, 0; -0.2, -0.3, 0];
f = 1000;
cpreal = [0.08,0.05,0.03;0.0408,0.03,0;0.03,0.02,0.03];
cpimag = [0,0.8,-0.004;-0.8,0,-0.01;0.004,0.01,0];
csound = 343.0;

E = TestEnergy(Xs,cpreal,cpimag,f,csound, micpos(:,1), micpos(:,2), micpos(:,3));
format long;
disp(E);
format short;