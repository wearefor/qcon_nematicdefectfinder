# qcon_nematicdefectfinder
MATLAB functions that identify topological half-integer defects in a nematic director field. Scripts use MATLAB’s conv2() function with a ring kernel to recursively perform line integrals that calculate the local topological index, producing a map that is then searched for regions of non-zero index. 

Files:
1. call_plotdefects.m  :  example usage
2. func_defectfind.m   :  core-function, identifies locations of half-integer defects provided a director field {nx, ny} in gridded form
3. func_defectorient.m :  identifies orientation of half-integer defects
4. func_gencircle.m    :  builds kernel for convolution
5. sampleframe         :  data used in example

Acknowledgements:
This research was conducted by Michael M. Norton (https://www.mmnorton.com) at Brandeis University, Waltham, Massachusetts USA in Seth Fraden's Lab (https://www.fradenlab.com/) where it was supported by the NSF MRSEC-1420382. Sample data was acquired by Linnea Lemma @ U.C. Santa Barbara (Zvonimir Dogic Lab http://dogiclab.physics.ucsb.edu/) and is featured in the work: Zhou et al., "Machine learning forecasting of active nematics," Soft Matter 2020 (https://pubs.rsc.org/en/content/articlelanding/2020/sm/d0sm01316a/)

