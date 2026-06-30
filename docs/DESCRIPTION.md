# Description of use cases

## Benchmark 1

This benchmark uses the data from [Oertlin et al.](https://academic.oup.com/nar/article/47/12/e70/5423604) (Supplementary File 1) that contains ground-truth for classification of regulatory layers.

Each tool is installed in its own environment and the *calling script* must accept the following arguments

```bash
run_tool.R --ribo RIBO --rna RNA --samples SAMPLES --out OUT --alpha ALPHA
```

Only these arguments are allowed. The values of `--ribo` `--rna`, and `--samples` are specified in the configuration file by `ribo`, `rna`, and `samples`, respectively, under `data`. The value of `--alpha` is the FDR threshold, given as `fdr_threshold` in the configuration file under `params`.

The output must be a csv table with the following columns

| Column name  | Description |
| ------------- | ------------- |
| identifier  | simulated gene ID  |
| class  | predicted class  |
| effect | log2 fold change or equivalent effect size |
| score | -log10(adjusted p-value) |

The class terminology is fixed according to the **anota2seq** labels. Hence, *translation* corresponds to *exclusive*, *abundance* to *forwarded*, and *buffering* to *buffered*. The *intensified* class is considered a sub-case of the mRNA abundance class, where the interaction term is non-zero and amplifies the RNA change. Features classified as *intensified* by **Ribotools** and **deltaTE** are thus re-labeled as *abundance*.

The choice of effect size and score depends on the methodology and class, *e.g.* for **Ribotools**, we use the interaction fold change and adjusted p-value, except for features in the *forwarded* and *buffered* classes, where the RNA fold change and adjusted p-value are used. For **anota2seq**, we use the native classification, given by `singleRegMode` (`anota2seqGetOutput` with `output="singleDf"`), to assign corresponding values for `apvEff` and `apvRvmPAdj`. Riborex provides no classification; all significant features are assigned to the *translation* class. Non significant features are assigned to the *background* class.

### Notes

* We use defaults for all tools, except otherwise stated.
* We use `fdr_threshold=0.1` and a default effect size of zero because neither **Riborex** nor the original implementation of the **deltaTE** method accept these parameters as input. DESeq2 is the default **Riborex** engine and is also used by **deltaTE**. In DESeq2, the default is to test that the log2 fold changes are equal to zero at a significance cutoff of 0.1.
* Effect size filtering `minEff` in **anota2seq** is applied post-hoc (after evaluating statistical significance), and is not directly comparable to the DESeq2 `lfcThreshold`. We set it to `minEff=NULL`.
* **anota2seq** and **Riborex** provide functions that must be integrated into an R script. **deltaTE** is not a tool *per se* (it cannot be installed); the script was taken from GitHub (master branch, commit 60b61ca, Apr. 2021) and patched to handle different conditions. **Ribotools** is a CLI-tool invoked with command-line options; the associated R script formats the input/output as expected, and calls **Ribotools** using `system2`. Values for `--alpha` and `--lfcThreshold` are passed via command-line arguments.
* The tested contrasts `control` *vs.* `treatment` use the default ordering of factor levels.

### Evaluation

Four test cases are defined: *(1)* Translation strict, *(2)* Translation only, *(3)* Any regulated, and *(4)* Buffering. Consult [run_evaluation.R](../workflow/scripts/run_evaluation.R) for details. This script is called with the following arguments

```bash
run_evaluation.R --rna RNA --translation TRAN --abundance ABUN --buffering BUFF --results FILE1 FILE2 [FILE3, ...] --out OUT
```

where `--rna` is the input count matrix, `--translation`, `--abundance`, and `--buffering` are the genes truly associated with these classes, given in the configuration file under `groundtruth`, and `--results` is a list of csv tables with expected columns (see above), one per tool tested.

We test *(A)* whether a tool's continuous score discriminates true positives from background, regardless of class assignment, and *(B)* the accuracy of a tool's explicit categorical classification.

The evaluation script should run if you add or remove tools, but keep the following in mind

* The input count matrices must have identical row names (gene or identifier).
* To ensure a fair comparison between tools that may differ in internal filtering, all tools are evaluated against the full set of input genes. Genes filtered out or not reported are assigned a score of zero and an effect size of zero, conservatively testing silent filtering as non-significant. For this data, this has no practical consequence, though, as all genes are expressed in all samples.
* Tool names are extracted from the paths given in `--results`. Any *fancy names* and their usage is hard coded.
* Tools that do not output classes for classification accuracy must be excluded for part *B* (hard coded).
* The figure layouts may have to be adjusted for a different number of tools (hard coded).

