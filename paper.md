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
  - name: Felipe Bordeu
    affiliation: 1
  - name: Julien Cortial
    affiliation: 1
  - name: Fabien Casenave
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

See Figure \autoref{fig:DeepLearningPrepost}.
![Example of deep learning prepost.\label{fig:DeepLearningPrepost}](DeepLearningPrepost.png){ width=20% }

# Acknowledgements

We acknowledge contributions from ...

# References
