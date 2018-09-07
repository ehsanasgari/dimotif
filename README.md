# DiMotif: Alignment-free Discriminative Protein Motif Miner

We present DiMotif as an alignment-free discriminative motif miner and evaluate the method for finding protein motifs in different settings. The significant motifs extracted could reliably detect the integrins, integrin-binding, and biofilm formation-related proteins on a reserved set of sequences with high F1 scores. In addition, DiMotif could detect experimentally verified motifs related to nuclear localization signals.


### DiMotif paper is currently under review and available on BioArXiv:
```
    @article {Asgari345843,
    author = {Asgari, Ehsaneddin and McHardy, Alice and Mofrad, Mohammad R. K.},
    title = {Probabilistic variable-length segmentation of protein sequences for discriminative motif mining (DiMotif) and sequence embedding (ProtVecX)},
    year = {2018},
    doi = {10.1101/345843},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2018/07/12/345843},
    eprint = {https://www.biorxiv.org/content/early/2018/07/12/345843.full.pdf},
    journal = {bioRxiv}
    }
```


<h1> User Manual </h1>

```
python3 dimotif.py --pos seqfile_of_positive_class --neg seqfile_of_negative_class --outdir output_directory --topn top_N_motifs --segs number_of_segmentations
```

Using the above mentioned command all the steps will be done sequentially and output will be organized in output directory.

<h3> Main parameters</h3>

--pos sequences file of the positive_class in txt or fasta format<br/>
--neg sequences file of the negative_class in txt or fasta format<br/>
--outdir output_directory <br/>
--topn how many motif to extract<br/>
--segs number of segmentation schemes to be sampled<br/>

## Computational Workflow

<ul>
<li>For a given set of positive sequences it extracts the most discriminative motifs in the positive class using a probabilistic segmentation inferred from Swiss-Prot</li>
<li>Motifs are hierarchically clustered according to their co-occurrence patterns in the positive sequences
Motifs are colored according to their most frequent secondary structure in PDB database</li>
<li>For each motif the normalized biophysical scores are also provided for further biophysical interpretations</li>
<li>The orange databases in the diagram are general-purpose databases and information. However, the red and blue databases are problem-specific datasets we want to find their related motifs.</li>
<ul>

![dimotif-2](https://user-images.githubusercontent.com/8551117/45029857-4d8db080-b04a-11e8-8e43-e42399c88217.png)


<img src='https://llp.berkeley.edu/wp-content/uploads/2018/07/Screen-Shot-2018-07-21-at-6.07.10-AM-1024x499.png'/>


