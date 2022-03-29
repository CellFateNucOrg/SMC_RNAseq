# SMC_RNAseq - ribo0

Code for processing ribo0 total RNA libraries from Moushumi Das with strains PMW382 (dpy-26::TEVcs;hs::TEV), PMW784 (scc-1::TEVcs;hs::TEV) and PMW775 (kle-2::TEVcs;hs::TEV) vs PMW366 (hs::TEV).

The general setup is similar to that described in the main branch, with a few differences:

- Salmon may not be the best aligner for repeats, so we also use STAR, BWA with HTSEQ to produce count tables.

- We used only simple contrasts of each strain to the TEV control (PMW366)

- Only one replicate was sequenced, so there is only laneNumber as a technical variable

- Repeat data needs to be downloaded from dfam.

- Repeats can be analysed as individual loci or as aggregates by repeat family - for the latter, use the scripts with rptAgg suffix
