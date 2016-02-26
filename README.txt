Description:
 This is a software package of the GMM-based compressive sensing (CS) inversion algorithm. The demo show the example of video CS. However, with minor modification, it can be used for image inpainting, image denoise, hyperspectral image CS.

Reference:
 Jianbo Yang Xuejun Liao, Xin Yuan, Patrick Llull, David J. Brady, Guillermo Sapiro and Lawrence Carin 
 ¡°Compressive Sensing by Learning a Gaussian Mixture Model from Measurements¡±  
 IEEE Transactions on Image Processing  (TIP)  vol. 24, no. 1, pp. 106-119, 2015.

Author:
 Jianbo Yang,  the department of Electrical and Computer Engineering, Duke University
 https://sites.google.com/site/jbysite/

Version history:
 version 1.0 (July 10, 2014)
 version 1.1 (Mar. 24, 2015)

Note:
1. To speed up algorithms MMLE-GMM and MMLE-MFA, the command "parafor" is used . If multi-core CPU is not available for your PC, you may want to change the files "fGMM_FR.m" and "fGMM_LR.m" accordingly.
2. To speed up algorithms MMLE-GMM and MMLE-MFA, we do not monitor the likelihood in the EM algorithm but simply set a maximum number of iterations in the parameter "para.iteMax". You may want to change "GMM_FR.m" and "GMM_LR.m" to monitor the likelihood in the EM algorithm.

Please send comments, suggestions and bug reports to email: jbyang.email@gmail.com


