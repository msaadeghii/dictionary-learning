%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  This package contains MATLAB functions and one demo implementing   %%%
%%%  the dictionary learning algorithms proposed in [1] and [2].                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Mostafa Sadeghi                                                                                                  %%%
%%%     Electrical Engineering Department,                                                                 %%%
%%%     Sharif University of Technology,                                                                       %%%
%%%     Tehran, IRAN.                                                                                                        %%%
%%%     E-mail: m.saadeghii@gmail.com                                                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a) Installation:

Prior to use this package, you need first to install the OMP and KSVD toolboxes
which are available at: http://www.cs.technion.ac.il/~ronrubin/software.html

Then, to install the current package, simply run "setup.m".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b) Package description:

PAU-DL.m        parallel atom-update dictionary learning algorithm proposed in [1] and [2]
OS-DL.m           one-stage dictionary learning algorithm proposed in [1]
APrU-DL.m      atom-profile update dictionary learning algorithm proposed in [1]
demo.m           a demo of dictionary recovery

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c) References:

[1] Mostafa Sadeghi and Massoud Babaie-Zadeh, and Christian Jutten,
    “Learning Overcomplete Dictionaries Based on Atom-by-Atom Updating,” 
    IEEE Trans. Signal Processing, vol. 62, n. 4, pp. 883-891, Feb. 2014.

[2] M. Sadeghi, M. Babaie-Zadeh, and C. Jutten, 
      “Learning over-complete dictionaries based on parallel atom-updating,” 
     in Proceedings of 23rd IEEE International Workshop on Machine Learning for Signal Processing (MLSP 2013), London, 2013.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%