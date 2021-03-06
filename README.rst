==============================
The python A-DaGO-Fun package
==============================
.. See the PDF package documentation for more information 
   on the use of the tool and different GO semantic simi-
   larity measures 
A-DaGO-Fun is a repository of python modules for analyzing 
protein or gene sets at the functional level based on Gene 
Ontology annotations using information content-based sema-
ntic similaritymeasures. It contains six main functions 
and implements 101 different functional similarity measures.

The main use cases of the library are:
  * Computing Information Content (IC), term and protein 
    semantic similarity scores: getTermFeatures, termsim 
    and funcsim.
  * Identifying enriched GO terms accounting for uncert-
    ainty in an annotation dataset: gossfeat.
  * Discovering functionally related or similar genes/pro-
    teins based on their GO terms: proteinfct.
  * Retrieving genes or proteins by their GO annotations 
    for disease gene and target discovery: proteinfit. 

Installation
------------

To install A-DaGO-Fun is quite straightforward and is similar 
to installation of any other python package. The whole package 
is relatively large (around 96Mb) and contains five modules 
and sets of files (GO term features and IC scores) available 
for download. 
 
Four packages, namely scipy, matplotlib, networkx and cPickle, 
need to be installed prior to the installation and use of 
A-DaGO-Fun. To install A-DaGO-Fun, the user needs to download 
the 'tar.gz' file and extract all files as follows::
::
tar -xzf dagofun.tar.gz

and then install from the `package' menu. To do this, one uses 
the following command::
::
python setup.py install --user

The A-DaGO-Fun package is free to use under GNU General Public 
License. You are free to copy, distribute and display information
contained herein, provided that it is done with appropriate ci-
tation of the tool. Thus, by using the A-DaGO-Fun package, it 
is assumed that you have read and accepted the agreement provi-
ded and that you agreed to be bound to all terms and conditions 
of this agreement. Please, use the following command line to 
see the package licence::
::
python setup.py --licence

The two commands above should be executed inside the dagofun
directory.
It is worth mentioning that A-DaGO-Fun is a portable python pa-
ckage and can be used without installing the whole package. You 
only need to be on a directory containing dagofun folder, which 
is a python library of A-DaGO-Fun.

Note that one module, namely tabulate.py for Pretty-print tabular 
data, which is borrowed from other authors, specifically written 
by `Sergey Astanin (s.astanin@gmail.com)' and collaborators.

Build status
------------

The main website for the A-DaGO-Fun package is 
http://web.cbio.uct.ac.za/ITGOM/adagofun where users can find 
essential information about A-DaGO-Fun. It is freely download-
able under GNU General Public License (GPL), pre-compiled for 
Linux version and protected by copyright laws. Users are free 
to copy, modify, merge, publish, distribute and display informa-
mation contained in the package, provided that it is done with 
appropriate citation of the package and by including the permi-
ssion notice in all copies or substantial portions of the modu-
le contained in this package.

It is currently maintained by one member of the core-develo-
pment team, Gaston K. Mazandu <gmazandu@gmail.com, 
gmazandu@cbio.uct.ac.za, kuzamunu@aims.ac.za, who regularly 
updates the information available in this package and makes 
every effort to ensure the quality of this information.

Administration
--------------

To start with A-DaGO-Fun package, type following commands:

  >>> import dagofun
  >>> help(dagofun)

Currently, this package provides six main modules: 
  * TermFeatures.py, TermSimilarity.py, ProteinSimilarity.py, 
  * ProteinIdentification.py, 
  * EnrichmentAnalysis.py  and 
  * ProteinClustering.py 
written independently and each containing a specific 
functions: 
  * getTermFeatures, termsim and funcsim, 
  * proteinfit, 
  * gossfeat and proteinfct, respectively. 

One can start a module GivenModule.py from the A-DaGO-Fun 
package as follows:

  >>> import dagofun.GivenModule as gm
  >>> help(gm)

After starting the module GivenModule.py as above, one can 
call or use the special function named gofunc of this module 
by writting

  >>> gm.gofunc(arguments)

The function named gofunc from the module GivenModule.py can 
also be made available directly as follows:

  >>> from dagofun.GivenModule import gofunc

After importing the function gofunc, to get help on how to use
the function, type the command above:

  >>> help(gofunc)

The function can be called directy as follows:

  >>> gofunc(parameters)

Finally, all the special functions in the package can be made 
available as follows:

  >>> from dagofun import *

and thus, each of the specific functions can be called directly 
as described above. Function arguments depend on the function 
and outputs a nicely formatted plain-text table, displayed either 
on the screen or written in a file, named according to the input 
data or file provided and the name displayed on the screen.

Detailed description and use of main functions
----------------------------------------------

(1) getTermFeatures
===================
This function from TermFeatures.py retrieves Information Content 
(IC) scores and other GO term features from the GO directed acyclic 
graph (DAG) structure. Given a GO term or list/tuple or a file 
containing a list of GO IDs, this function retrieves these charac-
teristics of these GO terms in the GO DAG, including their IC scores, 
provided the IC model. It uses GO-universal model ('u') by default, 
i.e., if no model is provided the GO-universal model is used.

(1.1) Usage:
------------
This function requires at least one argument and takes up to four
arguments.  
(a) a GO term, list/tuple or file containing a list of GO terms
(b) model, (c) drop and (d) output as described below.  

  >>> getTermFeatures(InputData, model='u', drop = 0, output=1)

This indicates that function argument (a) above is compulsory. 

* IC model:
Symbol of four different IC models implemeted in this package are:
   'u': For the GO-Universal (default if no IC model is provided)
   'w': For Wang et al.
   'z': For Zhang et al
   'a': For Annotation-based
See PDF documentation for more details about these IC models

* drop : boolean variable only useful in the context of Annotation-
based approach and it is set to 0 if Inferred from Electronic Anno-
tation (IEA) evidence code should be considered and to 1 otherwise.
 By default, it is set to 0.
     
*output: a boolean variable also set to 1 to display results on the 
screen and to 0 in a file. By default (output=1), results are disp-
layed on the screen, and finally default table display uses the pa- 
ckage module written by 'Sergey Astanin (s.astanin@gmail.com)' 
and collaborators. It used tablefmt="rst".
If results are written onto a file, the name of the file is basica-
lly the name of the first parameter in the function followed by TF
and  where ':' is replaced by '_', this is a case when this parame-
ter is a GO term.
	     
(1.2) Examples:
---------------
Note that you can be in any folder if the package was locally ins-
talled, otherwise make sure on a directory containing dagofun 
folder, which is a python library of A-DaGO-Fun.

  >>> from dagofun.TermFeatures import *
  >>> getTermFeatures('tests/TestTerms.txt', 'w', output=0)
  >>> terms = ['GO:0000001','GO:0048308', 'GO:0005385']
  >>> getTermFeatures(terms, 'a', 1)
  >>> getTermFeatures(terms, 'a', 1, 0)
  >>> getTermFeatures(terms, 'a', output=1)
  >>> getTermFeatures('tests/TestTerms.txt')
  >>> getTermFeatures('tests/TestTerms.txt', 'z')

Please check the file TestTerms.txt from the tests folder in the
dagofun directory. This file provides a model file of GO IDs.

(2) termsim
===========
This function from TermSimilarity.py computes semantic similarity 
scores between GO term pairs. Given a two GO IDs or a list/tuple
or two lists/tuples of GO IDs or a file containing a list of GO 
ID pairs, this function computes the semantic similarity scores
between GO ID pairs produced from these input GO IDs:

* For a given two GO terms as arguments, the function computes
  the semantic similarity between these two terms.
* For a list of tuple of GO IDs, the function computes similari-
  ty scores between all GO IDs pairs (a, b) for a and b in a 
  list or tuple with a != b.
* If two lists A and B are given, then similarity scores are 
  computed between all pairs (ai, bi) with ai in A and bi in B
  and 0 <= i <= min(len(A), len(B)) - 1
* If the file of GO ID pairs is provided, the function computes
  similarity scores between these GO ID pairs
   
(1.1) Usage:
------------
This function requires at least one argument and takes up to six
arguments.  
(a) A GO ID, list/tuple of GO IDs or file containing GO ID pairs
(b) A GO ID or a list/tuple or the ontology under consideration
(c) Approach to be used, the function supports up to for
    different approaches provided in a tuple/list. If non approach
    is provided, GO-universal model ('u') is used by default, 
(d) drop and (e) output as described previously.  

  >>> termsim('GOID1', 'GOID2', ontology='BP', approach='u', drop=0, output=1)
  >>> termsim(GO1, GO2, ontology='BP', approach='u', drop=0, output=1)
  >>> termsim(GO1, ontology ='BP', approach='u', drop=0, output=1)
  >>> termsim('FileName', ontology='BP', approach='u', drop=0, output=1)

* ontology: note that we are dealing with three ontology indepen-
  dently. Biological Process ('BP'), Molecular Function ('MF') and
  Cellular Component (CC). If no ontology is given then 'BP' is used

* Term Semantic Similarity (SS) Approach:
Approach or list/tuple of approaches under consideration (up to four 
(4) approaches can be considered). If more than 4 are given, these
beyond the 4th valid symbols are ignored.
   'u': for the GO-Universal
   'w': for Wang et al.
   'z': for Zhang et al
   'r': for Resnik
   'xr': for XGraSM-Resnik
   'n': for Nunivers
   'xn': for XGraSM-Nunivers
   'l':  for Lin
   'xl': for XGraSM-Lin
   'li': for Lin enhancement by Li et al (SimIC).
   's': for Lin enhancement by Relevance (Schlicker et al.)
See PDF documentation for more details about these term SS approaches
	     
(1.2) Examples:
---------------
Note that you can be in any folder if the package was locally ins-
talled, otherwise make sure on a directory containing dagofun 
folder, which is a python library of A-DaGO-Fun.

  >>> from dagofun.TermSimilarity import *
  >>> termsim('tests/TermSimTest.txt', approach = ('n','z','li'))
  >>> termsim('tests/TermSimTest.txt', output=0)
  >>> termsim('GO:0000001','GO:0048308', 'BP', ('xr', 'w', 'z', 'li'))
  >>> terms=['GO:0000001','GO:0048308', 'GO:0048311', 'GO:0000002']
  >>> termsim(terms, approach = ('n','z')) 

Please check the file TestSimTest.txt from the tests folder in the
dagofun directory. This file provides a model file of GO ID pairs.

(3) funcsim
===========
This function from ProteinSimilarity.py retrieves functional simi-
larity scores between protein/gene pairs. Given a string represen-
ting the name of the file containing protein IDs and their asso-
ciated GO IDs or two sets of GO IDs or a dictionary with protein
or gene IDs and keys and list ot tuples of GO IDs as values, for which 
functional similarity (SS) scores must be computed, this function 
computes the functional similarity scores between them:

* For a given file or an dictionary, if the target list of protein
or gene ID pairs is not provided, this function computes similari-
  ty scores between all protein/gene pairs (p, q) for p and q in  
  the file or dictionary.
   
(1.1) Usage:
------------
This function requires at least one argument and takes up to six
arguments.  
(a) Dictionary of protein/gene and  or file containing protein/
    genes and associated GO IDs annotating them.
(b) Protein or gene ID or a list/tuple or the ontology under consideration
(c) Approach to be used, the function supports up to for
    different approaches provided in a tuple/list. If non approach
    is provided, GO-universal model ('u') is used by default, 
(d) drop and (e) output as described previously.  

  >>> funcsim(GO1, GO2, ontology='BP', measure='ubma', drop=0, output=1)
  >>> funcsim(ProtGO, TargetPairs=[], ontology='BP', measure='u', drop=0, output=1)
  >>> funcsim(FileName, TargetPairs = [], ontology='BP', measure='ubma', drop=0, output=1)

* ontology: note that we are dealing with three ontology indepen-
  dently. Biological Process ('BP'), Molecular Function ('MF') and
  Cellular Component (CC). If no ontology is given then 'BP' is used

* FS measure
------------
measure or tuple of measures under consideration (up to three measures 
can be considered). The symbole of a given functional similarity 
measure is constructed as follows:
The starting letter r, n, l, li, s, x, a, z, w, and u represent GO term
 semantic similarity approaches and stand for Resnik, Nunivers, Lin, 
Li, Relevance, XGraSM, Annotation-based, Zhang, Wang and GO-universal, 
respectively. The suffixes gic, uic, dic, cou, cot, avg, bma, abm, bmm
hdf, vhdf and max represent SimGIC, SimUIC, SimDIC, SimCOU, SimCOT, 
Average, Best Match Average, Average Best Matches, Hausdorff, Variant
Hausdorff and Max measures, respectively. In cases where the
 prefix x is used, indicating XGraSM-based, it is immediately followed 
by the approach prefix. For example:
    'xlmax' for XGraSM-Lin based Average Functional Similarity Measure
    'agic' for Annotation-based SimGIC Functional Similarity Measure
    'zuic' for Zhang et al. based SimUIC Functional Similarity Measure
And the Union-Intersection functional similarity measure is represented
 by the symbol 'ui'.

- If the functional measure is not provided, the GO-universal based Best
Match Average (ubma) is used by default. 

See the PDF package documentation for more details about these FS 
measures and their corresponding symbols.
	     
(1.2) Examples:
---------------
Note that you can be in any folder if the package was locally ins-
talled, otherwise make sure on a directory containing dagofun 
folder, which is a python library of A-DaGO-Fun.

  >>> from dagofun.ProteinSimilarity import *
  >>> funcsim('tests/TestProteins.txt', measure = ('wcou','acou','acot'))
  >>> funcsim('tests/TestProteins.txt', ontology = 'BP', measure = 'agic', drop = 1, output=1)
  >>> funcsim('tests/TestProteins.txt', measure = ('wcou','acou','acot'), output=0)
  >>> proteinterms = {'Q5H9L2':['GO:0006355','GO:0006351'], 'P03891':['GO:0022904','GO:0044281','GO:0044237','GO:0006120'], 'Q5H9L2':['GO:0006355','GO:0006351']}
  >>> funcsim(proteinterms, measure = ('wvhdf','zvhdf','nvhdf'))
  >>> A = ['GO:0022904', 'GO:0044281', 'GO:0044237', 'GO:0006120']; B = ['GO:0006355', 'GO:0006351']
  >>> funcsim(A, B, measure = ('ub','nto','db','ub'))

Please check the file TestProteins.txt from the tests folder in the
dagofun directory. This file provides a model file of Protein-GO 
annotations.

(3) proteinfit
==============
This function from ProteinSearch.py retrieves genes or proteins con-
tributing to a given processes at a certain threshold or agreement 
level based of protein or gene annotations. For a given two strings 
representing the name of the file of background proteins, each with 
its GO ID annotations, and the target GO ID file or list, or a dic-
tionary with protein or gene IDs and keys and list ot tuples of GO 
IDs as values and the target GO ID file or list, the function iden-
tifies of candidate genes or proteins matching these terms or 
meeting the threshold or agreement level.

(3.1) Usage:
------------
This function requires at least two arguments and takes up to six
arguments.
(a) The name of the file containing the list of protein or gene IDs 
and their associated GO ID pairs. 
(b) A list/tuple or name of the file containing the list of GO ID 
targets for which, protein/gene should be identified.
(c) One of the GO ontologies: BP, MF and CC
(d) One of the term semantic similarity approaches. Refer to termsim
function
(e) The threshold score providing the semantic similarity degree at 
which terms are considered to be semantically close in the GO struc-
ture. 
(f) drop as described previously.
 
  >>> proteinfit(AnnotationData, TargetGOIDs, ontology='BP', approach='u', score=0.3, drop = 0)

The function outputs:
  (a) Summary statistics for different target GO IDs proteins, which 
includes following fields:
GO ID, GO term, Term Level, Number of proteins, Average SS, p-value 
and Bonferroni correction displayed on the screen or in a file 
depending on the argument output.
  (b) For each GO ID target, a summary statistics is provided a file 
named using the GO ID under consideration with following fields:
Protein ID, GO-ID related to the term,  GO-IDs with high-SS, 
Maximum SS and Average SS
  
A default parameters are GO-universal ('u') for SS, score = 0.3 the 
threshold or agreement level, Considering all GO evidence codes 
(drop = 0) and display by default on the screen (outputs=1), and 
finally default table display (tablefmt="rst") see tabulate package 
written by 'Sergey Astanin' and collaborators.  
	     
(3.2) Examples:
---------------
Note that you can be in any folder if the package was locally ins-
talled, otherwise make sure on a directory containing dagofun 
folder, which is a python library of A-DaGO-Fun.

  >>> from dagofun.ProteinSearch import *
  >>> targets = ['GO:0006355', 'GO:001905', 'GO:0001658']
  >>> proteinfit('tests/TestProteins.txt', target)
  >>> proteinfit('tests/TestProteins.txt', 'tests/TermSimTest.txt', 'BP', 'n')
  >>> background = {'Q5H9L2':['GO:0006355','GO:0006351'], 'P03891':['GO:0022904','GO:0044281','GO:0044237','GO:0006120'], 'Q5H9L2':['GO:0006355','GO:0006351']}
  >>> proteinfit(background, targets, approach='s', score = 0)

Please check the file TermSimTest.txt and TestSimTest.txt from the 
tests folder in the dagofun directory. This file provides a model 
file of GO ID pairs.

(4) gossfeat
============
This function from EnrichmentAnalysis.py biological retrieves 
processes most pertinent to the experiment performed based on the
 target set and background provided. Given two strings represen-
ting the name of the file of background proteins, each with its 
GO ID annotations, and the target protein file, the function 
identify biological processes most pertinent to the experiment 
performed.
The function incorporates the complex dependence structure of the 
GO DAG and the uncertainty in annotation data using fuzzy expres-
sions through GO term semantic similarity measures.

(4.1) Usage:
------------
This function requires at least two arguments and takes up to six
arguments.
(a) The name of the file containing the list of protein or gene IDs
and their associated GO ID pairs, as described previuosly.
(b) The name of the file containing the list of target proteins or 
genes.
(c) Ontology, approch, threshold score and drop arguments work in 
the same way as described previously (see 3.1, and 2.1).
(d) The significance level cut-off (pvalue) from which an identified 
term is considered to be statistically significant, set to 0.05 by 
default.

  >>> gossfeat(ReferenceFile, TargetFile, ontology='BP', approach='u', score=0.3, pvalue=0.05, drop=0)

It is worth mentioning that the strict cases where the score 0 and 
1 correspond to the well known traditional cases: score=0 when 
considering true path rule and score=1 for effective occurence or
exact match

(4.2) Examples:
---------------
  >>> from dagofun.EnrichmentAnalysis import *
  >>> gossfeat('tests/ReferenceSetTest.txt', "tests/TargetSetTest.txt", approach = 's', score=0.5)
  >>> gossfeat('tests/ReferenceSetTest.txt', "tests/TargetSetTest.txt")
  >>> gossfeat('tests/ReferenceSetTest.txt', "tests/TargetSetTest.txt", approach = 'n', score=0.0)
  >>> gossfeat('tests/ReferenceSetTest.txt', "tests/TargetSetTest.txt", approach = 'w', score=1.0)
Note that for this specific function, resulting list all enriched 
terms and their features are displayed in a file from the forder 
where the module is being executed. For the above example, the 
name of this file will be ReferenceSetTestEA.txt

Please check the file ReferenceSetTest.txt and TargetSetTest.txt from 
the tests folder in the dagofun directory. This file provides a model 
file of GO ID pairs.

(5) proteinfct
==============
This function from ProteinClustering.py allows the partitioning of a 
gene or protein set into a set of biological meaningful sub-classes 
using their functional closeness based on GO annotations and derived 
from a selected semantic similarity model. For a given two strings 
representing the name of the file of proteins, each with its GO ID 
annotations, and the potential target proteins file or list to be 
clustered, or a dictionary with protein or gene IDs and keys and list 
ot tuples of GO IDs as values and the target protein file or list, the 
function elucidates functionally related or similar genes/proteins 
based on their GO termsidentifies of candidate genes or proteins 
matching these terms or meeting the threshold or agreement level.

(5.1) Usage:
------------
This function requires at least two arguments and takes up to nine
arguments.
(a) Dictionary of protein/gene and  or name of the file containing 
protein/genes and associated GO IDs annotating them.
(b) Protein or gene ID or a list/tuple or the ontology under conside-
ration
(c) For the functional similarity measure, refer to funcsim (see 3.1)
(d) The threshold score providing the functional similarity degree at 
which proteins can be considered to be functionally close, set to 0.3
by defaut. The score 0 indicates that all protein pairs with func-
tional similarity score greater that 0.
(d) drop and output as described previously.  
(e) Clustering model under consideration, and this function implements
three different models (mclust):
  - hierarchical clustering (mclust = 1)
  - Graph spectral clustering or kmeans (mclust = 2)
  - community detecting model by Thomas Aynaud, 2009 (mclust = 3)
(f) Number of clusters (nclust) applies only for the kmeans model and
it is set to 0 by default. In this case, if mclust is less than 2 then
the comminity detecting model is applied instead of kmeans!

  >>> proteinfct(AnnotationData, TargetIDs=[], ontology='BP', measure='ubma', score=0.3, mclust=1, nclust=0, drop=0, output=1)

(5.2) Examples:
---------------
  >>> from dagofun.ProteinClustering import *
  >>> proteinfct('tests/TestProteins.txt', measure='nbma', score=0.0, mclust=2, nclust=3)
  >>> proteinfct('tests/TestProteins.txt', mclust=2)
  >>> proteinfct('tests/TestProteins.txt', measure='rbmm', score=0.0, mclust=1)
  >>> background = {'Q5H9L2':['GO:0006355','GO:0006351'], 'P03891':['GO:0022904','GO:0044281','GO:0044237','GO:0006120'], 'Q5H9L2':['GO:0006355','GO:0006351']}
  >>> proteinfct(background, measure='rbmm', score=0.0)

Version history
---------------

- 15.1: Initial A-DaGO-Fun release in June 2015.

Package URL
-----------

http://web.cbio.uct.ac.za/ITGOM/adagafun

Maintainer
----------

Gaston K. Mazandu
Email: gmazandu@cbio.uct.ac.za, gmazandu@gmail.com, 
       kuzamunu@aims.ac.za

Contributors
------------

Gaston K Mazandu, Emile R Chimusa, Mbiyavanga Mamana, Nicola J Mulder
Emails: gmazandu@cbio.uct.ac.za, emile@cbio.uct.ac.za, 
        mamana@aims.ac.za, nicola.mulder@uct.ac.za
