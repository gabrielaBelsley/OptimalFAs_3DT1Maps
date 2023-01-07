# Optimal FAs for 3D T1 mapping using VIBE

The scripts implement the theoretical framework described in the paper:

**Optimal Flip Angles for In Vivo Liver 3D T1 Mapping and B1+ Mapping at 3T**

Gabriela Belsley (1), Damian J. Tyler (1), Matthew D. Robson (1,2), Elizabeth M. Tunnicliffe (1)

(1) Oxford Centre for Clinical Magnetic Resonance Research, Division of
Cardiovascular Medicine, Radcliffe Department of Medicine, University of Oxford,
Oxford, UK.

(2) Perspectum, Oxford, UK

---------------------------------------------

Run first the script varianceT1SPGR_BeforeScanInVivo.m.

This script calculates the T1 Coefficient of Variation (CoV) for a range of T1s, B1+ factors, VIBE SNRs, B1+ SNRs and different FA combinations.


Then run script varianceT1SPGR_OptimalFAs.m that finds the optimal FAs using a min-max approach. 
The simulation figures from the paper (Figures 1-3) and Table 1 are generated with this script.
Note: also download the folder DrosteEffect-BrewerMap-ca40391 in order to generate the figures.

This code is distributed under the MIT license. If you use it, please cite the code: 

Author: Gabriela Belsley, University of Oxford, gabriela.belsley@stx.ox.ac.uk
OR gabi.belsley@gmail.com

Please contact me, if you have any questions. 

---------------------------------------------

