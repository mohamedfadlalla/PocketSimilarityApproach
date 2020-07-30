# Pocket Similarity Approach
This repository help reporduce the results of the manuscript at:
https://chemrxiv.org/articles/preprint/COVID19_Approved_Drug_Repurposing_Pocket_Similarity_Approach/12722483

this script does not include fpocket and pocketmatch algorithms; you need to do this from the orginal software repository:
Pocketmatch 2.1:http://proline.physics.iisc.ernet.in/pocketmatch/
Fpocket: https://github.com/Discngine/fpocket

the script works fine on python 3.7.1

dependancies:
1. Biopython
2. Pandas
3. Numpy

Drug_pocket_prepration.py:
handle extracting the pocket of the drug

COVID19_pocket_prepration.py:
handle extracting the pocket COVID19

