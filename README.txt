(1) GenEv_refit.C takes Output_GenEv.root to produce:
Output_GenEv_1.root - needed variables,
Refit.root - contains fitted values, calculated missing mass, probability and chi2.
(2) macro dohist.C executes hist.C (1 command instead of 3); it takes Refit.root in order to produce hist.root (various histograms) and some png's.