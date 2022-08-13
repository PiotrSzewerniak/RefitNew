Each measurement, regardless of the quantity we reconstruct, is biased with a certain error in determining this value. It is obvious due to the lack of a  device with perfect accuracy and efficiency.

The main task of the kinematic refit procedure is the best possible transformation of our measurements (using known uncertainties) in such a way that their values are as close as possible to the expected values of these quantities. For this purpose, two mathematical methods are combined: the Least Squares Method and the Method of Lagrange Multipliers. The first one allows one to calculate the minimal differences between the value obtained as a result of the measurement and the value calculated by us using the uncertainties.  The second one allows us to impose certain limitations on the measurements (like missing or invariant  mass condition), which we can choose ourselves. These additional conditions are defined by the conservation laws of physics.

Obtained results, used data and all "how-to" are presented in my Master Thesis: https://www.overleaf.com/read/phqytysxmhtx

Code would not have been created without help from D.Sc. Izabela Ciepał and PhD Bogusław Włoch.

In order to launch program you need to have ROOT Cern installed.
(1) GenEv_refit.C takes Output_GenEv.root to produce:
Output_GenEv_1.root - needed variables,
Refit.root - contains fitted values, calculated missing mass, probability and chi2.
(2) macro dohist.C executes hist.C (1 command instead of 3); it takes Refit.root in order to produce hist.root (various histograms) and some png's.
