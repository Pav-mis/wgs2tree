# wgs2tree
This is a nextflow implementation of a workflow for the rapid generation of consensus phylogenomic trees from whole genome sequence data.

The workflow is as follows:

1. Samples, either consisting of a single assembly, or multiple assemblies (ie, hap1 and hap2) are used as input
2. For each assembly, single-copy orthologues are identified with BUSCO
3. For samples consisting of more than one assembly, Buscomp is used to compile the maximal set of Busco orthologues
4. For each retrieved BUSCO gene, a multiple sequence alignment is performed with MAFFT
5. A phylogenetic tree is generated for each gene with IQTree2
6. A consensus phylogenomic tree is generated with ASTRAL
