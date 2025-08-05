# Refinement-expansion substitution model (resub)

Resub is a class of amino acid substituion models that assumes there were 19 amino acids at the top of the tree and 20 at the bottom. This accounts for expansions and refinements to the coding alphabet in early evolutionary history. 



## Installation instructions


This package requires BEAST 2.7.8. or newer.

1. Launch BEAUti
2. Click on `File` -> `Manage Packages`
3. Install `resub`. 


### Other dependencies


Requires BEAGLE installed to use the fast multi-epoch substitution models. 

To run these XML files, please also install the following BEAST 2 packages through the package manager
- BEASTLabs
- Gamma spike model



## Running

This model currently does not have BEAUti support, however there are several `XML` files in the `examples` directory that can help you get started. These files can also be used to reproduce the results of the article.


## Support

BEAST user forums [https://groups.google.com/g/beast-users](https://groups.google.com/g/beast-users)

Jordan Douglas jordan.douglas@auckland.ac.nz


## References

Douglas, J., Bouckaert, R., Carter Jr, C. W., & Wills, P. R. (2025). Reduced amino acid substitution matrices find traces of ancient coding alphabets in modern day proteins. Molecular Biology and Evolution (in press). Preprint: https://www.biorxiv.org/content/10.1101/2025.02.19.639152v1 







