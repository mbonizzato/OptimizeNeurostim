# OptimizeNeurostim

GP-BO optimization of neurostimulation.

Companion code for the following paper : [Link to paper available upon publication]

Please cite the present code as:
Bonizzato M., Guay-Hottin R., Laferrière S., Caron-Grenier O., Lajoie G. and Dancause N. 2021. “OptimizeNeurostim”, GitHub.  https://github.com/mbonizzato/OptimizeNeurostim/.


### Data

The dataset can be downloaded at :  [Link to data available upon publication]

Please cite as:
Bonizzato M., Massai E.&, Côté S.&, Quessy S., Martinez M., and Dancause N. 2021. “OptimizeNeurostim” OSF. osf.io/54vhx.

&, these authors equally contributed to the dataset.



## Custom Matlab implementation

Tested with all Matlab version between 2018b and 2021b.
Does not require specialized third-party libraries.


## GPML Matlab implementation

Tested with all Matlab version between 2018b and 2021b.
Requires GPML Matlab Code version 4.2 by Carl Edward Rasmussen and Hannes Nickisch
http://www.gaussianprocess.org/gpml/code/matlab/doc/

The software included the following line, where "evalc" is used to suppress screen output.

    evalc('hyp = minimize(hyp, @gp, -10, infprior, [], covf, likf, x, y);');

For the sake of performance, we recommend substituting this line with:

    hyp = minimize(hyp, @gp, -10, infprior, [], covf, likf, x, y);

and manually modifying the GPML minimize function to suppress all screen outputs.


## Python implementation


To install Gpy (https://sheffieldml.github.io/GPy/):

    pip install gpy
   
Requires a recent version of scipy (1.3.0 or later) and numpy.


We suggest replacing the linalg.py file located in GPy\util with the one in this repository to avoid issues with the jitchol method (https://github.com/SheffieldML/GPy/issues/660). We used np.linalg.cholesky instead of scipy.linalg.cholesky and raised the number of regularization attempts. We also removed a warning in force_F_ordered about the arrays not being F order.

We noticed that about 5% of the value of hyperparameters with uniform prior (kernel lengthscales and variance) were outside of the prior's support. This behavior might be related to https://github.com/SheffieldML/GPy/issues/639.

In the python version, we used a bounded constraint on the model's Gaussian noise hyperparameter instead of a smooth uniform prior as in the matlab implementation to avoid numerical instability.


