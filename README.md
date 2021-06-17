# Blind-deconvolution-based-on-the-ratio-of-cyclic-content
This is the Matlab code of the blind deconvolution based on the ratio of cyclic content (BD-RCC). BD-RCC can be used to recover repetitive impacts signal in noisy vibration.

The simulation example shows that BD-RCC superior than CYCBD, MCKD, MOMEDA and MCG-L1/L2-D.

MCG-L1/L2-D：Minimum correlated generalized L1/L2 deconvolution     
MCKD：Maximum correlated kurtosis deconvolution      
CYCBD：Maximum second-order cyclo-stationarity blind deconvolution      
MOMEDA：Multipoint optimal minimum entropy deconvolution adjusted       

The main difference between various deconvolution methods lies in the characterization index.There are several indicators with similar definitions to RCC, and the details are as follows:

Some equivalent definitions of RCC can be found in papers [1-3].

Hilbert envelope spectrum fault feature ratio (HESFFR) [1]

Envelop harmonic-to-noise ratio (EHNR)[2-3]

Note：This procedure is only for manuscript review. For other purposes, please contact aresmiki@163.com

[1]	W. He, Y. Zi, B. Chen, F. Wu, Z. He, Automatic fault feature extraction of mechanical anomaly on induction motor bearing using ensemble super-wavelet transform, Mech. Syst. Signal Process. 54–55 (2015) 457–480. https://doi.org/10.1016/j.ymssp.2014.09.007.

[2]	M. Zhao, X. Jia, A novel strategy for signal denoising using reweighted SVD and its applications to weak fault feature enhancement of rotating machinery, Mech. Syst. Signal Process. 94 (2017) 129–147. https://doi.org/10.1016/j.ymssp.2017.02.036.

[3]	X. Xu, M. Zhao, J. Lin, Y. Lei, Envelope harmonic-to-noise ratio for periodic impulses detection and its application to bearing diagnosis, Measurement. 91 (2016) 385–397. https://doi.org/10.1016/j.measurement.2016.05.073.


