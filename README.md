# Estimation of Sparse Functional Additive Models with Adaptive Group LASSO

The folder consists of simulations and analysis for both tecator data and air pollutant data.

1.File - smoothing.R: file consists of codes for two functions

GCV.fit: implement generalized cross-validation for smoothing spline,

est.smooth: employ smoothing spline for estimation.

File - simulation ssFAM.R: This file consists of codes for our simulation studies.

Folder - "tecator analysis"

3.a: tecator.txt: the tecator data, its description is available on http://lib.stat.cmu.edu/datasets/tecator; 

3.b: matlab code for PACE: files used to implment pace for both training data and test data;

3.c: traintec.csv and testtec.csv: they contain the results from the matlab code used to implment pace for both training data and test data;

3.d: Tecator analysis.R: data analysis for tecator data. Different estimation methods are compared.

Folder - "pm 2.5 analysis"

4.a: pm25_2000.csv: the functional predictors in the training data and test data, which is obtained from data cleaning.R,

4.b: matlab code for pace: files used to implment pace for both training data and test data;

4.c: train.mat and test.mat: they contain the results from the matlab code used to implment pace for both training data and test data;

4.d: airpollution.R: data analysis for air pollutant data. Different estimation methods are compared.
