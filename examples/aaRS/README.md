

# Instructions


This directory contains 2 aaRS datasets. The model is configured in `config.json` file, which can be injected into the `xml` file using the `-df` option. 

- `resub_aaRS.xml` contains the common elements of the catalytic domain (Class I and II)
- `resub_uryzme.xml` contains a subset of the sites in the alignment above (Class I and II)


The two Class trees are linked onto the same timescale using divergence time estimates from a previous tree of life study. 

1. To select the cherry, edit `config.json` to contain the right cherry, and also the right MRCA node that the alphabet transition should occur at. For example to use the `IV` cherry make sure that `state1=I` and `state2=V`, and that the `focalPrior` points to the common ancestor of IleRS and ValRS:

```
{
 "state1":"I",
 "state2":"V",
 "focalPrior":"IleValRS.prior"
}
```

The other focal priors are: 
`AlaGlyRS.prior` for the `AG` cherry
`AspLys.prior` for the `DK` cherry
`AspAsnLysRS.prior` for the `DN` cherry
`GlnGluGlxRS.prior` for the `EQ` cherry
`PheHisRS.prior` for the `FH` cherry
`ThrProRS.prior` for the `PT` cherry
`SerGlyRS.prior` for the `SG` cherry
`TrpTyrRS.prior` for the `WY` cherry


2. To run resub:
```
~/beast/bin/beast -df config.json resub_aaRS.xml
```

This assumes that BEAST 2 is installed in the `~/beast/bin/` directory, and that all the right packages are installed, plus BEAGLE.

For BEAGLE GPU acceleration, use
```
~/beast/bin/beast -df config.json -beagle_gpu resub_aaRS.xml
```



