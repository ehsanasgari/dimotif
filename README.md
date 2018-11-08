# DiMotif: Alignment-free Discriminative Protein Motif Discovery

We present DiMotif as an alignment-free discriminative motif discovery method and evaluate the method for finding protein motifs in three different settings: (1) comparison of DiMotif with two existing approaches on 20 distinct motif discovery problems which are experimentally verified, (2) classification-based approach for the motifs extracted for integrins, integrin-binding proteins, and biofilm formation, and (3) in sequence pattern searching for nuclear localization signal. The DiMotif, in general, obtained high recall scores, while having a comparable F1 score with other methods in discovery of experimentally verified motifs. Having high recall suggests that the DiMotif can be used for short-list creation for further experimental investigations on motifs.


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

<h1> DiMotif step-by-step </h1>

An ipython notebook containing an example of motif discovery using DiMotif is provided here: 
https://github.com/ehsanasgari/dimotif/blob/master/notebook/DiMotif_step_by_step_example.ipynb



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


