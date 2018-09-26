clear all;

load('TestData.mat');

vlb = [-1 -1 0 0 -1 -1 0 0];
vub = [+1 +1 2 0.5 +1 +1 2 0.5];

XA = mic_info(:,1).';
YA = mic_info(:,2).';
ZA = mic_info(:,3).';

% Crossover rate
pc = 0.75;
% Multiplication factor
F = 0.4;
% Number of runs
Nruns = 1;
q = 128;
% Number of generations
Ng = 250;

cpimag = imag(CSM);
cpreal = real(CSM);
f = frequencies;

tic;
for k3 = 1 : Nruns
    fprintf('\tf = %d Hz, run %d\n', f, k3);
    seed_nr = abs( sum(100*clock) - cputime );
    rand('state',seed_nr) % set the random generator to a different state each time
    clear seed_nr

    % run DE
    pop_total = [];
    energy_total = [];
    E1 = zeros(1, q);
    % Population size
    startpop = lhd_new(q, vlb, vub);
    X1 = startpop;
    % energy function values for elements of X1
    for k4 = 1:q
        E1(k4) = Energy_fun_C(X1(k4,1:end),cpreal,cpimag,f,X1(k4,end),XA,YA,ZA);
    end
    % tranfer to column vector
    E1 = E1';
    pop_total = [pop_total; X1];
    energy_total = [energy_total E1'];
    for k5 = 1:Ng-1
        % create_descendants_DE_1
        % OUTPUT:
        % dummy1 = popinterm, i.e. the partner population
        % n = size of partner population
        % X3 = popnew, i.e., the new population obtained from
        % crossover between popold and popinterm
        [dummy1, Nt, X3] = create_descendants_DE_1(X1, F, pc, vlb, ...
                                                   vub, pop_total);
        % new_generation_DE_1
        % OUTPUT:
        % X4 = popnew, i.e. the population for the next generation
        % dummy2 are the corresponding values for the energy
        % function
        % pop_total and energy_total are updated with respect to
        % input values
        [X4, dummy2, pop_total, energy_total] = ...
            new_generation_DE_1('Energy_fun_C', cpreal, cpimag, ...
                                f, c, XA, YA, ZA, X1, E1, ...
                                X3, vlb, vub, pop_total, energy_total);
        E1 = dummy2;
        % update initial population
        X1 = X4;
        % save convergence behaviour
        [G2(k3,k5), Imin] = min(E1);
        F2(:,k3,k5) = X1(Imin,:)';
    end
    [~, Imin] = min(E1);
    F1(:,k3) = X1(Imin,:)';
    G(k3) = E1(Imin);

end
toc