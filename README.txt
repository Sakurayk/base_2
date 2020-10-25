Sample MATLAB codes: Fast singular value shrinkage using CPA.
(c) Masaki Onuki, 2017

Author:
Masaki Onuki (masaki.o@msp-lab.org)

These MATLAB codes are the examples of fast singular value shrinkage using Chebyshev polynomial approximation (CPA). When you use this samples for your paper, please cite the paper below:
M. Onuki, S. Ono, K. Shirai, and Y. Tanaka, ``Fast Singular Value Shrinkage with Chebyshev Polynomial Approximation Based on Signal Sparsity,'' IEEE Transactions on Signal Processing.
This paper includes theoretical and practical details of the fast singular value shrinkage with CPA.

Prerequisites:
MATLAB 2015b or later. We have not confirmed yet if our codes can be run by using the older versions than MATLAB 2015b.

Instruction:
We prepared the application indicated at Section V-F in the above paper. In the ``Sample_code_Fast_SVS_CPA'', the singular value shrinkage using CPA and the exact method was applied to the matrix rank minimization. You can select the CPA-based method or the exact method by changing the variable ``CPA_SVS''. The CPA-based method is used if CPA_SVS = 1, or the exact method is used if CPA_SVS = 0. The matrix rank and the missing rate of used data can be determined by changing the variables ``rank'' and ``Miss_rate''. In addition, choose the approximation order for CPA by changing the variable ``Approx order''. Finally, run the source code. Enjoy!

Update history:
Aug. 17, 2017: v0.1 - (original release)