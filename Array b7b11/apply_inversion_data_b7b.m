close all;
clear all;
clc;
fprintf('*** Benchmark case b7b-Sarradj optimization ***\n\tReading data...\n');
CsmReaderb7b;
num_freq = length(freqs);

% Now run DE for finding the source position that provides the CSM with a
% maximum match between model predictions ans provided CSM
% source position X,Y,Z, and amplitude.^2
vlb = [-0.25 -0.25 0 0 -0.25 -0.25 0 0 -0.25 -0.25 0 0 -0.25 -0.25 0 0 csound];
vub = [0.25 0.25 1.5 2 0.25 0.25 1.5 2 0.25 0.25 1.5 2 0.25 0.25 1.5 2 csound];

XA = mic_info(1,:);
YA = mic_info(2,:);
ZA = mic_info(3,:);

% Crossover rate
pc = 0.75;
% Multiplication factor
F = 0.4;
% Number of runs
Nruns = 50;
q = 128;
% Number of generations
Ng = 600;

fprintf('\tStart optimzation for %d frequencies...\n', num_freq);
for fi = 61:62
    f = freqs(fi);
    
    tic;
    fprintf('\tFreq. %d Hz: ', f);
    reverseStr = '';
    for k3 = 1 : Nruns
        msg = sprintf('Run %d/%d...', k3, Nruns);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
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
            E1(k4) = Energy_fun_C(X1(k4,1:end-1),cpreal,cpimag,fi,f,csound,XA,YA,ZA);
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
                                    fi , f, csound, XA, YA, ZA, X1, E1, ...
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
    msg = sprintf('Runs completed! Total time: %d s\n', round(toc));
    fprintf([reverseStr, msg]);
    
    fprintf('\tSaving for f = %d Hz\n\n', f);
    nowstr = datestr(now,'HH.MM_dd-mm');
    save(['.\DATA\b7b_' num2str(f) 'HZ.mat'], ...
        'G', 'G2', 'F1', 'F2', 'pc', 'F', 'Nruns', 'q', 'Ng', 'f', ...
        'vlb', 'vub');
    
end
fprintf('*** Finished ***\n');