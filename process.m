%Import our original data and divide it into sub block. For each block, rebuild it by recover.m

%{
    input£ºk£ºAssumed sparsity level of the original data. Used to determine pool numbers.
           pcell£ºApproximately equal to the number of cells used in each pool.
           p£ºDetermining the columns of M together with pcell.
           alpha£ºNoise. When the value is 1, there is no noise. 
    output£ºX£ºThe original data of each sub block
            M£ºThe measurement matrix for each sub block
            R£ºThe inference of each sub block
%}

function [X,M,R] = process(k,pcell,p,alpha)

original = importdata('count4070.txt');
XX = [original.data];
XX = XX';

ncolM = floor(pcell/p);	% the number of columns per M

C = {};

% If the column size of the last block is bigger than half of the previous one, 
% it will become a new seperate block, else merge it into previous one.
ni = floor(size(XX,1)/ncolM);
for i = 1:ni
    C{i} = XX((ncolM*(i-1)+1):ncolM*i,:);
end
if size(XX,1)-ncolM*ni>ncolM/2
    C{ni+1} = XX(ncolM*ni+1:size(XX,1),:);
else
    C{ni} = [C{ni};XX(ncolM*ni+1:size(XX,1),:)];
end


M = {};% Used to save the measurement matrix of each block
R = {};% Used to save the inference of each block
X = {};% Used to save the original data of each sub block
for i = 1:length(C)
    XX = C{i};
    X = [X XX];
    [MM,recoverXX] = recover(floor(k*size(XX,1))*2,size(XX,1),p,XX,alpha);
    M = [M MM];
    R = [R recoverXX];
end

