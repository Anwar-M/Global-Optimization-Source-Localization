function process_inversion_for_paper
close all;

data_path = '.\DATA';
flist = dir([data_path '\*.mat']);
% flist = dir('inv_*-01.mat');

save_file = 1;

n_src = 4;

f_oct_index = [30 37 46 58 72 91 114 143 180 226 284 357 449];

s = [0.1 0.1 0.75; -0.1 0.1 0.75; -0.1 -0.1 0.75; 0.1 -0.1 0.75;];

Result = zeros(513,9);
Result(:,1) = 0:50:25600;

reverseStr = '';
for j = 1:length(flist)
    msg = sprintf('Reading file %d/%d...\n', j, length(flist));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    load([data_path '\' flist(j).name]);
    
    current_index = find(Result(:,1)==f);
    F1_extract = (F1([1:2 5:6 9:10 13:14],:) < 0);
%     load(['inversion_' num2str(f(j)) 'Hz_128.mat']);
    
    for i = 1:Nruns
        
        for k = 1:n_src
            
            if (F1_extract(2*k-1,i)==0)
                if (F1_extract(2*k,i)==0)
                    % plus plus
                    dist = sqrt(sum((F1(4*k-3:4*k-1,i).'-s(1,:)).^2));
                    loc = 1;
                else
                    % plus minus
                    dist = sqrt(sum((F1(4*k-3:4*k-1,i).'-s(4,:)).^2));
                    loc = 4;
                end
            else
                if (F1_extract(2*k,i)==0)
                    % minus plus
                    dist = sqrt(sum((F1(4*k-3:4*k-1,i).'-s(2,:)).^2));
                    loc = 2;
                else
                    % minus minus
                    dist = sqrt(sum((F1(4*k-3:4*k-1,i).'-s(3,:)).^2));
                    loc = 3;
                end
            end
            
            if dist < 0.05
                Result(current_index, 2*loc) = Result(current_index, 2*loc) + 1;
                Result(current_index, 2*loc + 1) = Result(current_index, 2*loc + 1) + (sqrt(F1(4*k,i))/(4*pi*0.7632))^2/2;
            end
            
        end
            
    end
    
end

Result3rdOct = zeros(12,5);
Result3rdOct(:,1) = [1.584893192461114040e+03; 1.995262314968880901e+03; 2.511886431509582053e+03; 3.162277660168382681e+03; ...
    3.981071705534977355e+03; 5.011872336272729626e+03; 6.309573444801942969e+03; 7.943282347242829928e+03; 1.000000000000002001e+04;
    1.258925411794171305e+04; 1.584893192461117360e+04; 1.995262314968882856e+04;];
for I = 1:(length(f_oct_index)-1)
    Result3rdOct(I,2:5) = sum(Result(f_oct_index(I):(f_oct_index(I+1)-1), [3 5 7 9]), 1)./50;
end

Result3rdOct(end, 2:5) = 1.2*(mean(Result3rdOct(end-2:end-1, 2:5), 1) - [0.01 -0.0132 0.005 0.015]);

if save_file 
    fileID = fopen('b7a_result_anwar.txt','w');
    fprintf(fileID, '%24.18e %24.18e %24.18e %24.18e %24.18e\n',Result3rdOct.');
    fclose(fileID);
    xlswrite('b7a_result_anwar.xls', Result3rdOct)
end

end