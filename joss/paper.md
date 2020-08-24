---
title: 'Ripserer.jl: flexible and efficient persistent homology computation in Julia'
tags:
  - Julia
  - Persistent Homology
  - Topological Data Analysis
authors:
  - name: Matija Čufar
    affiliation: 1
affiliations:
 - name: Independent Researcher
   index: 1
date: 24 August 2020
bibliography: paper.bib
---

# Summary

Persistent homology [@edelsbrunner2008persistent] is a relatively recent computational
technique that extracts topological information from various kinds of datasets. This
topological information gives us a good overview of the global shape of the data as well as
giving us a description of its local geometry. Since its introduction it has been used in a
diverse range of applications, including biology [@bernoff2016biological], material
science [@lee2017quantifying] to signal processing [@tralie2016high] and computer
vision [@asaad2017topological]. A problem persistent homology faces is the very large size of
combinatorial structures it has to work with. Recent algorithmic advances employ various
computational shortcuts to overcome this problem.

Among the most successful implementations of persistent homology is
Ripser [@bauer2019ripser]. With its speed and low memory usage, it makes persistent homology
practical for larger datasets, even in higher dimensions. The introduction of Ripser has
spawned a whole cottage industry of extensions and wrappers. Some examples include
Ripser++ [@zhang2020gpu], Lock-free Ripser [@morozov2020towards], Ripser.py [@tralie2018ripser]
and Cubical Ripser [@kaji2020cubical].

Julia [@bezanson2017julia] is a relatively recent language specifically designed for
technical computing. Even though the language has only been stable since 2018, it already
offers packages that are among the best in various areas such as differential
equations [@rackauckas2017differentialequations], mathematical optimization [@dunning2017jump]
and network analysis [@bromberger2017lightgraphs]. The main advantages it offers are a high
degree of composability through multiple dispatch, and performance that can compete with
statically compiled languages like C, C++ and Fortran.

In this work, we present Ripserer.jl, a pure Julia implementation of persistent homology
based on the algorithm that powers Ripser. The main goal of this project is to provide a
highly generic and extensible persistent homology library that is easy to use for
non-experts without sacrificing too much performance.

# Package Details

Along with its companion package PersistenceDiagrams.jl, Ripserer.jl provides a featureful
environment for computing persistent homology and integrating it with the rest of Julia's
data science stack. At the time of writing, it offers the following features.

* Fast Vietoris-Rips, alpha complex and cubical persistent homology computation.
* Representative cocycle and critical simplex computation.
* Support for coefficients in any, possibly user-defined, field.
* Convenient persistence diagram and representative cocycle visualization via
  Plots.jl[^1] recipes.
* Bottleneck and Wasserstein matching and distance computation.
* Various persistence diagram vectorization functions implemented with persistence
  images [@adams2017persistence] and persistence curves [@chung2019persistence].
* Easy extensibility through a documented API.

![Example visualizations. The plot on the left shows the three main representative cocycles
in the data. The right plot show the persistence diagram.](figure.png)

Ripserer.jl provides a way to hook into its persistent homology algorithm through an API at
various stages in the algorithm. This allows researchers to experiment with different
approaches to persistent homology without having to reimplement the algorithm from scratch
or forking an existing repository. The extensions can range from simple sparsification
schemes to complicated reimplementations, like the included cubical homology extension.

Being written in the high-level programming language Julia is also an advantage as
understanding the codebase does not require knowing C++, which is a notoriously difficult
language to learn.

Our benchmarks[^2] show that Ripserer's performance is very close to that of Ripser. It
tends to be slightly slower for dense inputs and slightly faster for very sparse inputs. In
the cubical case, we compared it to Cubical Ripser. There the performance was worse, taking
up to 3 times as long to compute some results. This is expected as Cubical Ripser is much
more specialized for its use case and even splits its code into different repositories for
different dimensions.

We have not compared performance with newer, parallel implementations such as Ripser++ or
lock-free Ripser. Judging from the benchmarks they provide, we expect them to perform much
better. Their downside, however, is that they require powerful hardware, such as GPUs or
large numbers of processors.

[^1]: https://github.com/JuliaPlots/Plots.jl
[^2]: https://mtsch.github.io/Ripserer.jl/dev/benchmarks/

# Acknowledgments

We would like to thank Žiga Virk for comments, suggestions and ideas, and Ulrich Bauer for
making the source code of Ripser freely available.

# References
