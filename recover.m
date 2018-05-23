%For each sub block, we generate a measurement matrix M and add noise e. 
%Using the noisy measurement matrix eM together observed value y to inference recoverXX.

%{
    input£ºm£ºThe number of rows for each sub block.
           n£ºThe number of columns of each sub block.
           p£ºDetermine the proportion of 1 in the M matrix.
           XX£ºOriginal data. Is used to generate the observed value y.
           alpha£ºNoise. When the value is 1, there is no noise.
    output£ºM£ºthe measurement matrix of each sub block.
            recoverXX£ºinference sub X.
%}

function [M,recoverXX] = recover(m,n,p,XX,alpha)
recoverXX = [];
M = double(rand(m,n)<p);  %Use uniform distribution of [0,1] to randomly generate M (m¡Án).
                          %For each element in M, if value is less than p, we set to 1, otherwise 0.
                                                                     

%The code below is used to generate the full observation matrix.
%{ 
a(1:round(n/4))=1;
a(round(n/4)+1:round(2*n/4))=0.1;
a(round(2*n/4)+1:round(3*n/4))=0.01;
a(round(3*n/4)+1:n)=0.001;
M=[];
for i=1:m
    a = a(randperm(length(a)));
    M = [M;a];
end
%}

e = unifrnd(alpha,2-alpha,m,n); %Use uniform distribution of [alpha,2-alpha] to randomly generate erorr(m¡Án).
eM = e.*M;                                  
for i = 1:size(XX,2)
	y = eM*XX(:,i);            
	if sum(y) == 0              
		xp = 0*XX(:,i);         
		recoverXX = [recoverXX xp];
	else
		[w, hcond] = linsolve(M*M',y);
		x0 = M'*w;
		xp = l1eq_pd(x0, M, [], y);
		recoverXX = [recoverXX xp];
	end
end
end

