---
title: 'BasicTools: a numerical simulation toolbox'
tags:
  - Python
  - C++
  - mesh
  - fields
  - finite elements
  - post-treatment
authors:
  - name: Felipe Bordeu^[first author]
    affiliation: 1
  - name: Julien Cortial
    affiliation: 1
  - name: Fabien Casenave^[corresponding author]
    orcid: 0000-0002-8810-9128
    affiliation: 1
affiliations:
 - name: Safran Tech, Digital Sciences & Technologies Department, Rue des Jeunes Bois, Châteaufort, 78114 Magny-Les-Hameaux, France
   index: 1
date: 2 September 2021
bibliography: paper.bib

---

# Summary

Numerical simulations of physical phenomena are computed by various available codes, 
but many similar operations are required when producing these simulations. Post-treatments
of physical fields are a common need, as well as handling and modifying meshes and 
fields. BasicTools is a library addressing these supporting tasks. It contains an 
efficient data model for mesh and field objects and input/output routines compatible 
with various formats. A finite element engine allows the assembling of abstract 
variational formulation, differential operators and the integration of fields on 
volume and surfaces.

# Statement of need

BasicTools has been used in various projects in artificial intelligence and 
model order reduction for physical problems [@ROM-net], [@mca26010017], 
[@UQindustrialDesign], [@datatargetVAE], topological optimization (travaux 
citant BasicTools ?) and material sciences (travaux BigMeca).



# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Figures

![Example of deep learning prepost.\label{fig:DeepLearningPrepost}](DeepLearningPrepost.png){ width=100% }
See Figure \autoref{fig:DeepLearningPrepost}.


# Acknowledgements

We acknowledge contributions from ...

# References

<!-- 
# Citations
Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->