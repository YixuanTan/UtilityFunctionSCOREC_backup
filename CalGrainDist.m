function [Gnum, Gsize] = CalGrainDist(dimension, lowbound, upbound, steps)
stepstr = num2str(steps, '%04d');
file = strcat('voronmc.', stepstr,'.txt');
A = importdata(file);
condition1 = ((A(:,1) >= lowbound ) | (A(:,1) < upbound ));
A = A(condition1, :); % keep the rows satisfying condition. 
B = A(:,3);
C = sortrows(B);
D = [];
index = C(1);
Didx = 1;
D(1,1) = 0;

for k = 1:size(C,1)
	if index == C(k)
	  D(Didx) = D(Didx) + 1;
	else
	  index = C(k);
	  Didx = Didx + 1;
	  D(Didx) = D(Didx) + 1;
	end
	D(Didx+1) = 0;
end
	
GrainSize = (D(1:end-1))';
size(GrainSize)
sum(GrainSize);
mean(GrainSize)
