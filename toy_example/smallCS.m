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
	x0 = M'*y; 
	xp = l1eq_pd(x0, M, [], y);
	recoverX = [recoverX xp];
end

fname = ['pool_',num2str(pool),'_p_',num2str(p),'.mat'];
%save fname -ascii recoverX
save(fname,'recoverX','M','-mat');

