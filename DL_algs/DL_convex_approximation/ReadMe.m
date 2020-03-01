%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   MATLAB implementations of the convex relaxation approach for dictionary
%   learning proposed in the following paper:
%
%  Mostafa Sadeghi, Massoud Babaie-Zadeh, and Christian Jutten, 
%  “Dictionary learning for sparse representation: A novel approach,” 
%  IEEE Signal Processing Letters, vol. 20, no. 12, pp. 1195-1198, December 2013.
%
%  version 1, Feb. 26, 2014.
%  
%  Prior to use, please download and install the OMP and K-SVD toolboxes
%  from https://www.cs.technion.ac.il/~elad/software
%  
%  IMPORTANT:
%  Please carefully read Section III. A (last sentences of the first paragraph)
%  of the above mentioned paper about modifying part of the algorithms in the
%  case of synthetic or real experiments. In other words, according to the
%  application at hand, you should (un-)comment:
%  line 433 in "NewMOD.m", 
%  line 447 in "NewMDU.m", and 
%  line 438 in "NewSGK.m"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Mostafa Sadeghi
%     Electrical Engineering Department,
%     Sharif University of Technology,
%     Tehran, IRAN.
%     E-mail: m.saadeghii@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
