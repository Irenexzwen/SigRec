# SigRec
Compressed sensing scheme recovering single cell expression profile 
https://doi.org/10.1101/338319

\\

Put ***star.m***, ***process.m***, ***recover.m***, ***l1eq_pd.m*** in the same folder. You only need use **MATLAB** to run ***star.m***.

Enter the following code in the command-line window.

```matlab
star(k,pcell,p,alpha);
```


***star.m*** calls ***process.m***.（line21）and calculates the correlation coefficient.

***process.m*** divides the data into block and calls ***recover.m*** to infer the data. (line42)

(**WARNING**) In ***process.m***, you need to set the path of your data file. (line15)

And You need to make sure that the read data matrix is <code>cells×genes</code>.

***recover.m*** calls ***l1eq_pd.m*** to complete Compressed Sensing.(line43)

You can use ***smallCS.m*** as an example to run.
