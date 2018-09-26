function [popnew, val, pop_total, energy_total] = ...
    new_generation_DE_1(fname, cpreal, cpimag, fi, f, csound, XA, YA, ZA, ...
                        popold, fvalold, ui, vlb, vub, pop_total, ...
                        energy_total)

NP = size(popold,1);    % number of population elements (q)
popnew = popold;
val = fvalold;

%------------------------------------------------------------------------------------------------------
% Select which vectors are allowed to enter the new population (compare popold and ui)
%------------------------------------------------------------------------------------------------------
for k=1:NP
  tempval = feval(fname,ui(k,1:end-1),cpreal,cpimag,fi,f,ui(k,end),XA,YA,ZA); % check cost of competitor  ui

  if (abs(tempval) <= abs(fvalold(k)))    % if competitor is better than value in "cost array"
     popnew(k,:) = ui(k,:);               % replace old vector with new one (for new iteration)
     val(k)   = tempval;                  % save value in "cost array"
  end
  
  % save all evaluated vectors with their energies (10 march 2004, Camiel)
  pop_total = [pop_total; ui(k,:)];
  energy_total = [energy_total tempval];
  clear tempval
end

% end of program