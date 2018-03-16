Thank you for using this package for Granger causality estimation.

This is a Matlab GUI for estimating Granger causatliy on time series data set. It is still under test, so any feedback from you would be helpful. Free for acadamic use only.

* version log by Qiang on 10th Sep 2015
The current version 2.0.0.1 is tested in Matlab2014a, and the previous version 1.0.0.1 was developed and tested in Matlab2008.  



EASY TO RUN

In the current working path, the GUI would become avaliable if you type 
>>GrangerCausalityGUI
in the command window of Matlab.
Then the use of this tool should be straightforward.


METHODS

We provide three different methods to estimate Granger causality, including Conditional Granger Causality, Partial Granger Causality, spatio-temporal Granger Causality, and Granger causality with signal-dependent noise.

In spatio-temporal Granger causality (stGC), we need to specify if we want to use 
   1) optimal time window division on the time series to estimate the averaged GC across time windows, 
or 2) equal time windows with fixed window numbers
or 3) no time window division.

We also need to specify the spatial resolution. For example, if we are calculating the GC between two ROIs with m and n voxles, respectively. We can specify the voxels for each ROI in the GUI and compute the GC between all pairs of voxels from two ROIs. Then, the GC between these two ROIs can be estimated by the mean GC among all paris of voxels. 

In Granger causality with signal-dependent noise (GCSDN), we need to specify the model order in variance, ie. BEKK order.
The frequency domain result is also avialable for GCSDN. The function for detecting the existence of SDN has also been included in this package, but it will work for time series data with repeat blocks only. 



RESULTS

The results are shown in tables and figures. The results are also saved in 'myresult.mat' in the current working path. 

For temporal GC,
causMatrix(j,k)   ----causality from k-th time series to j-th time series
pvalMatrix(j,k)   ----p-value of estimated  causality, if avaliable
TtimeCau{j,k}      ----causality from k-th time series to j-th time series in each time window

For spatial GC and spatio-temporal GC, 
causMatrix(j,k)   ----mean causality from k-th ROI to j-th ROI, which is averaged over all pairs of voxels from two ROIs
stdMatrix(j,k)    ----std in the estimated causality

For conditional and partial GC
meanArray(j,k)  ----mean causality from k-th time series to j-th time series given by bootstrap
varArray(j,k)   ----std causality from k-th time series to j-th time series given by bootstrap

For GC with signal-dependent noise
TcausMatrix(j,k) --- causality from k-th time series to j-th by GCSDN
TpvalMatrix(j,k) --- p-value of estimated causality, if avaliable

CONTACT

This package is built by Jianfeng Feng group, in accompany of the papers "Spatio-temporal Granger Causality: A New Framework. NeuroImage, 2013, 79:241-63.”, ”Uncovering interactions in the frequence domain. PLoS Comput Biol 2008, 4(5): e1000087”, and "Attention-dependent Modulation of Cortical Taste Circuits Revealed by Granger Causality with Signal-dependent Noise. PLoS Computational Biology, 2013, 9(10): e1003265"
Contact: mrqiangluo@gmail.com
Last modified by Qiang Luo at 12:22 2/28/2014.