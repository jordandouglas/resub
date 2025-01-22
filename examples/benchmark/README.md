
# Instructions


This directory contains 16 benchmark datasets used to screen the effects of resub. These datasets are each specified in a `json` file, which can be injected into the `resub.xml` file using the `-df` option. The datasets named `BB...` are taken from BAliBASE 3.0, and the remainders are aminoacyl-tRNA synthetases.


1. Navigate to a folder, e.g.,
```
cd BB50008
```

2. Select which pair of amino acids (i.e., cherry) to operate on by editing the `metadata.json` file. For example, to work on the `IV`  cherry, make sure that `state1=I` and `state2=V`:


```
{
 "dataset":"BB50008",
 "N": 26,
 "L": 810,
 "Lcore": 211,
 "desc":"Papain-like protease (pfam PF00112)",
 "taxonomy":"Eukaryota,Viruses",
 "state1":"I",
 "state2":"V",
 "data":"<sequence totalcount='20' ..."
}
```


3. To run resub on this dataset/cherry combination:
```
~/beast/bin/beast -df metadata.json ../resub.xml
```

This assumes that BEAST 2 is installed in the `~/beast/bin/` directory, and that all the right packages are installed, plus BEAGLE.

For BEAGLE GPU acceleration, use
```
~/beast/bin/beast -df metadata.json -beagle_gpu ../resub.xml
```



