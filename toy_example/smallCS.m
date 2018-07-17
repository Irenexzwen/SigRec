
%{
    input：pool：The number of pools.
           p：Determine the proportion of 1 in the M matrix.
    output：recoverXX：inference X.
%}

function smallCS(pool,p)
pool = str2num(pool)
p = str2num(p)

original = importdata('pan_expression.txt');
X = [original.data];

%Transform to random binary matrix  by value comparison logic computation
M = double(rand(pool,64)<p); 

recoverX = [];
for i = 1:size(X,2)
	y = M*X(:,i);          
	[w, hcond] = linsolve(M*M',y);
	x0 = M'*w; 
	xp = l1eq_pd(x0, M, [], y);
	recoverX = [recoverX xp];
end

fname = ['pool_',num2str(pool),'_p_',num2str(p),'.mat'];
save(fname,'recoverX','M','-mat');
