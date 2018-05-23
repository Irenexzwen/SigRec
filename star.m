%{
    input£º
		   k£ºAssumed sparsity level of the original data. Used to determine pool numbers.
           pcell£ºApproximately equal to the number of cells used in each pool.
           p£ºDetermining the columns of M together with pcell.
		   
		   
	out£º  A mat file include:	
			   X(The original data of each sub block)
			   M(The measurement matrix for each sub block)
			   R(The inference of each sub block)
			   median correlation of all genes 
%}

function star(k,pcell,p,alpha)
k = str2num(k);
pcell = str2num(pcell);
p = str2num(p);
alpha = str2num(alpha);

[X,M,R] = process(k,pcell,p,alpha);


%Correlation analysis
XX = [];
recoverXX = [];	% Combine each sub block to original dimension
for i = 1:length(R)
    XX = [XX;X{i}];
    recoverXX = [recoverXX;R{i}];
end

recoverXX(recoverXX<0) = 0;	% Cleaning our reconstructed data
recoverXX = round(recoverXX);

CC=[];
for i = 1:size(XX,2)
    C = corrcoef(XX(:,i),recoverXX(:,i));	
    C = C(1,2);
    CC = [CC;C];
end

data = CC;
data(isnan(data)) = 1;


an = mean(data);
dian = median(data);


fname = ['k_',num2str(k),'_pcell_',num2str(pcell),'_p_',num2str(p),'_a_',num2str(alpha),'_pj_',num2str(an),'_zw_',num2str(dian),'.mat'];
save(fname,'X','M','R','an','dian','-mat');
%end 
