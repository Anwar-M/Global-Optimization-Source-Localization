function [popinterm,n,popnew] = create_descendants_DE_1(popold, F, CR, vlb, vub, pop_total);

% popnew has as many vectors as popold; each vector of popnew lines up with the same vector
% popinterm has in general more vectors than popold

NP = size(popold,1);      % number of population elements (q)
D  = size(popold,2);      % number of variables
n  = 10*NP;                % number of vectors: taken at the save side, part of them will be thrown away due to doubles

% Matrix A consists of rows of 3 integers, indicating the indices of vectors in the old population
% Each integer has a value between 1 and q. The number of rows is n.
% The new vector is calculated as:   vector(index1) + F*vector(index2) -
% F*vector(index3)
% Create A (n rows, 3 columns). The factor 0.9999999999999 makes sure that
% matrix A contains no value of NP + 1
A = 1+floor(NP*0.9999999999999*rand(n,3));
% make sure that per row all elements are different
[I,J] = find(A(:,1)~=A(:,2) & A(:,1)~=A(:,3) & A(:,2)~=A(:,3));
% A can now have any length, < n
A = A(I,:); clear I J
% make sure that all rows are different. The first column is sorted (from
% low to high)
A = unique(A,'rows');
% update n
n = length(A(:,1));

% the first row of A can be e.g. 1 3 5. Matrix M then becomes: [1 0 F 0 -F zeros(1,q-5)].
% This matrix M is used for calculating the partners.
% ~ (1:n)' are the row indices
% ~ A(:,1), A(:,2), and A(:,3) are the column indices.
M = zeros(n,NP);
M(sub2ind(size(M),(1:n)',A(:,1))) = 1;
M(sub2ind(size(M),(1:n)',A(:,2))) = F;
M(sub2ind(size(M),(1:n)',A(:,3))) = -F;
% calculate all n partners: each row of popinterm contains a partner
% popinterm is of size (n,D)
popinterm = M*popold; % each row of popinterm contains a potential descendant 
clear A M

% select admissable vectors
% (crossover on admissable vectors does not alter admissability)
lowerbound = ones(n,1)*vlb;
upperbound = ones(n,1)*vub;
z1 = all(popinterm >= lowerbound,2);
z2 = all(popinterm <= upperbound,2);
z = z1 .* z2;
I = find(z==1);
if isempty(I)
    popinterm = popold;
else
    popinterm = popinterm(I,:);
end
clear lowerbound upperbound z1 z2 z I

n = length(popinterm(:,1));

for p = 1 : NP
    % apply crossover from 1 vector of popold with all vectors of popinterm
    mui = rand(n,D) < CR;         % all random numbers < CR are 1, 0 otherwise
    mpo = mui < 0.5;              % inverse mask to mui: mui added to mpo gives a matrix consisting of ones
    % popinterm2 does not result in other parameter VALUES than those in pop_old and popinterm, but in other
    % combinations of the D parameter values
    popinterm2  = (ones(n,1)*popold(p,:)).*mpo + popinterm.*mui;
    % popinterm2 = setdiff(popinterm2, pop_total, 'rows'); % Camiel, 10 march 2004
    % select at random one element from the admissable set
    n2 = length(popinterm2(:,1));
    if (n2 == 0)
        disp('line 79 in create_descendants_1.m')
        keyboard
    end
    element = randperm(n2);
    popnew(p,:) = popinterm2(element(1),:);
    clear mui mpo popinterm2 n2 element
    
end % for p = 1 : NP

% end of program