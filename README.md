# GP-EST: A Unifying Algorithm for GP-UCB and Probability of Improvement
This is a matlab demonstration for an algorithm for Bayesian optimization with the Gaussian process assumption. The algorithm is fully described in

Optimization as Estimation with Gaussian Processes in Bandit Settings (Zi Wang, Bolei Zhou, Stefanie Jegelka), In International Conference on Artificial Intelligence and Statistics (AISTATS), 2016.

This paper is available at http://zi-wang.com/pub/wang-aistats16.pdf.

To run the code, first install the gpml toolbox by Carl Rasmussen and Hannes Nickisch (http://www.gaussianprocess.org/gpml/code/matlab/doc/). See gpo_example.m for a full explanation on how to run the algorithms.

In the paper, we also conducted experiments related to trajectory optimization and image classification. The trajectory optimization experiment is based on the Airplane2D example the Drake toolbox at http://drake.mit.edu. The image classification experiment follows the paper "Learning deep features for scene recognition using places database (Zhou et. al.)" in NIPS 2014. And all the datasets are available online.




