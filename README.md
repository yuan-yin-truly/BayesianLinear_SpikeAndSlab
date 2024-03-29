[![DOI](https://zenodo.org/badge/342649594.svg)](https://zenodo.org/doi/10.5281/zenodo.10648564)
# A Bayesian Linear Model of Multiple High-throughput Sequencing Data under Unknown Environmental Conditions

### Abstract

The regulatory relationship between transcription factors (TF) and target gene (TG) can alter under different environmental conditions, and it is typically hidden to us what experiment conditions correspond to which regulatory relationship. Therefore, even with a large number of public RNASeq data labeled with experiment conditions, uncovering the regulatory relationships are not straightforward. ChIPSeq experiments inform us the TFs that bind near TG. For ChIPSeq setup that artificially induces the production of a TF and saturate cells with the TF, we may assume those TFs that bind near TG represents the entire regulator space.

If we treat each TF and TG as a random variable, we may summarize the data we have as follows: an unlabeled collection of instantiations of random variables from multiple unknown data generating processes, although we do know all the TFs whose (unknown) subset participated in each data generating process. We present a statistical model to tackle this problem: find clusters in RNASeq data, select true regulators from TFs in ChIPSeq data for each cluster using a spike- and-slab prior strategy, and run a linear regression of TG vs. true regulators in each cluster. These subtasks are built in a single model, and we develop a Metropolis-within-Gibbs sampling strategy for inference on the parameters of the model.

*For the full report, please refer to [__report.pdf__](https://github.com/yuan-yin-truly/BayesianLinear_SpikeAndSlab/blob/main/report.pdf) at repo __BayesianLinear\_SpikeAndSlab__.*
