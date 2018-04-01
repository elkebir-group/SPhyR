We consider simulated datasets and a real metastatic colorectal cancer dataset [[1]](#ref1).

The simulated datasets are obtained in three steps.

1. We use the `ms` program
[[2]](#ref2) to [generate](perfect_phylogeny/generate.sh) perfect phylogeny trees, which are in the directory [perfect_phylogeny](perfect_phylogeny).
2. We [generate](k_dollo/generate.sh) losses using the
`simulate` executable, the resulting matrices are in the
directory [k_dollo](k_dollo).
3. We [generate](flip/generate.sh) errors using the `perturb`
executable yielding the matrices in directory [flip](flip).

## References
<a name="ref1"></a>
[1] Leung, M. L., Davis, A., Gao, R., Casasent, A., Wang, Y., Sei, E., et al. (2017). Single cell DNA sequencing reveals a late-dissemination model in metastatic colorectal cancer. [*Genome Research*, gr.209973.116.](http://doi.org/10.1101/gr.209973.116)
<a name="ref2"></a>
[2] Hudson, R. R. (2002). Generating samples under a Wright–Fisher neutral model of genetic variation. [*Bioinformatics*, 18(2), 337–338.](http://doi.org/10.1093/bioinformatics/18.2.337)
