# SigRec
Compressed sensing scheme recovering single cell expression profile 
https://doi.org/10.1101/338319

Put ***star.m***, ***process.m***, ***recover.m***, ***l1eq_pd.m*** in the same floder. You only need run ***star.m***.

***star.m*** calls ***process.m*** 

~~~matlab
```[X,M,R] = process(k,pcell,p,alpha);```
~~~

and calculates the correlation coefficient.

***process.m*** divides the data into block and calls ***recover.m*** to infer the data. 

~~~matlab
```[MM,recoverXX] = recover(floor(k*size(XX,1))*2,size(XX,1),p,XX,alpha);```
~~~

(**WARNING**) In ***process.m***, you need to set the path of your data file. 

~~~matlab
```original = importdata('count4070.txt');```
~~~

And You need to make sure that the read data matrix is $cells×genes$.

​***recover.m*** calls ***l1eq_pd.m*** to complete Compressed Sensing.

~~~matlab
```xp = l1eq_pd(x0, M, [], y);```
~~~

You can use ***smallCS.m*** as an example to run
