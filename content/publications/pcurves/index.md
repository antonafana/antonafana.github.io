---
title: 'Principal Curves In Metric Spaces And The Space Of Probability Measures'

# Authors
authors:
  - Andrew Warren
  - me
  - Forest Kobayashi
  - Young-Heon Kim
  - Geoffrey Schiebinger

date: '2025-05-07'

# Schedule page publish date (NOT publication's date).
publishDate: '2017-01-01T00:00:00Z'

# Publication type.
# Accepts a single type but formatted as a YAML list (for Hugo requirements).
# Enter a publication type from the CSL standard.
publication_types: ['article']

# Publication name and optional abbreviated publication name.
publication: On arXiv
publication_short: On arXiv

abstract: We introduce principal curves in Wasserstein space, and in general compact metric spaces. Our motivation for the Wasserstein case comes from optimal-transport-based trajectory inference, where a developing population of cells traces out a curve in Wasserstein space. Our framework enables new experimental procedures for collecting high-density time-courses of developing populations of cells - time-points can be processed in parallel (making it easier to collect more time-points). However, then the time of collection is unknown, and must be recovered by solving a seriation problem (or one-dimensional manifold learning problem). We propose an estimator based on Wasserstein principal curves, and prove it is consistent for recovering a curve of probability measures in Wasserstein space from empirical samples. This consistency theorem is obtained via a series of results regarding principal curves in compact metric spaces. In particular, we establish the validity of certain numerical discretization schemes for principal curves, which is a new result even in the Euclidean setting. 

# Summary. An optional shortened abstract.
summary: We take a fresh look at the principal curves model in probability space and apply it to scRNA-seq time courses.

tags:
  - High-dimensional
  - scRNA-seq
  - curve-fitting

# Display this page in the Featured widget?
featured: true

# Standard identifiers for auto-linking
hugoblox:
  ids:
    doi: https://doi.org/10.48550/arXiv.2505.04168

# Custom links
links:
  - type: pdf
    url: "https://arxiv.org/pdf/2505.04168"
  - type: code
    url: "https://github.com/antonafana/wasserstein_principal_curves"
---
