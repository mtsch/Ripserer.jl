# Benchmarks

The following tables show benchmarks that compare Ripserer's performance with
[Ripser](https://github.com/Ripser/ripser), [Cubical
Ripser](https://github.com/CubicalRipser/), and
[Eirene.jl](https://github.com/Eetion/Eirene.jl). The benchmarking code and more info about
the datasets are available [here](https://github.com/mtsch/RipsererBenchmarks.jl).

All benchmarks were performed on a laptop with an Intel(R) Core(TM) i5-4200U CPU @ 1.60GHz
with 8GB of RAM.

We used [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl/) to perform the
timing benchmarks and [Valgrind's Massif
tool](https://www.valgrind.org/docs/manual/ms-manual.html) to measure peak heap sizes
(i.e. total memory footprint).

The benchmarks were performed with Ripserer v0.15, `master` versions of Ripser (commit
hash `286d369`) and Cubical Ripser (commit hashes `6edb9c5` for 2D and `a063dac` for 3D),
and Eirene v1.3.5.

The timings show the minimum time taken among five runs of the benchmark.

The heap sizes for Ripserer include the Julia runtime.

## Comparison with Ripser

In this experiment, we performed benchmarks with the datasets presented in the [Ripser
article](https://arxiv.org/abs/1908.02518). We only used the datasets that we were able to
run with less than 8GB memory. All datasets were parsed as `Float32` as that is what Ripser
supports. The time it takes to parse a file is included for both Ripser and Ripserer.

### Dense results

|dataset       |size|dim|threshold|Ripserer|Ripser   |ratio|Ripserer heap|Ripser heap|
|:-------------|:---|:--|:--------|:-------|:--------|:----|:------------|:----------|
|`o3_1024`     |1024|3  |1.8      |4.576 s |3.057 s  |1.497|374.1 MiB    |151.0 MiB  |
|`o3_4096`     |4096|3  |1.4      |151.527 s|76.177 s|1.989|4.7 GiB      |4.1 GiB    |
|`dragon2000`  |2000|1  |         |3.133 s |2.833 s  |1.106|316.7 MiB    |296.8 MiB  |
|`fract-r`     |512 |2  |         |22.807 s|19.482 s |1.171|2.2 GiB      |2.0 GiB    |
|`random16`    |50  |2  |         |8 ms    |10 ms    |0.803|111.1 MiB    |1.1 MiB    |
|`sphere_3_192`|192 |2  |         |1.549 s |1.491 s  |1.039|287.0 MiB    |209.5 MiB  |

### Sparse results

These benchmarks were performed with the `sparse=true` keyword argument.

|dataset       |size|dim|threshold|Ripserer|Ripser   |ratio|Ripserer heap|Ripser heap|
|:-------------|:---|:--|:--------|:-------|:--------|:----|:------------|:----------|
|`o3_1024`     |1024|3  |1.8      |3.036 s |3.057 s  |0.993|418.2 MiB    |151.0 MiB  |
|`o3_4096`     |4096|3  |1.4      |76.052 s|76.177 s |0.998|4.9 GiB      |4.1 GiB    |
|`dragon2000`  |2000|1  |         |3.588 s |2.833 s  |1.267|350.4 MiB    |296.8 MiB  |
|`fract-r`     |512 |2  |         |25.399 s|19.482 s |1.304|2.2 GiB      |2.0 GiB    |
|`random16`    |50  |2  |         |9 ms    |10 ms    |0.932|111.1 MiB    |1.1 MiB    |
|`sphere_3_192`|192 |2  |         |1.734 s |1.491 s  |1.163|288.5 MiB    |209.5 MiB  |

### Alpha-Rips

These benchmarks were performed on sparse matrices that correspond to the 1-skeleta of
Delaunay triangulations. The purpose of these is to show performance with very sparse
inputs.

|dataset              |size |dim|Ripserer  |Ripser    |ratio|Ripserer heap|Ripser heap|
|:--------------------|:----|:--|:---------|:---------|:----|:------------|:----------|
|`alpha_3_sphere_3000`|3000 |3  |636 ms    |789 ms    |0.807|138.4 MiB    |33.2 MiB   |
|`alpha_torus_10_000` |10000|2  |872 ms    |1.179 s   |0.741|130.0 MiB    |27.7 MiB   |
|`alpha_5_sphere_1000`|1000 |5  |49.431 s  |46.707 s  |1.058|387.2 MiB    |202.0 MiB  |
|`alpha_dragon_2000`  |2000 |2  |56 ms     |76 ms     |0.744|2.4 GiB      |1.5 GiB    |
|`alpha_4_sphere_2000`|2000 |4  |5.844 s   |6.203 s   |0.942|110.9 MiB    |33.2 MiB   |

## Comparison with Cubical Ripser

In these benchmarks, we used some of the datasets presented in the [Cubical
Ripser](https://arxiv.org/abs/2005.12692) article. We limited the 2D image size to 1999×999
as the current `master` (commit hash `6edb9c5`) version of 2D Cubical Ripser throws an
assertion error for anything larger. We were also unable to perform 3D 256×256×256 image
benchmarks due to Ripserer running out of memory. The `eltype` of all datasets is `Float64`,
because that is what Cubical Ripser supports. When running Ripserer in the real world, it's
a good idea to use the image's native data types. This will _slightly_ reduce the memory
footprint and increase performance.

|dataset       |size   |dim|Ripserer  |Cubical Ripser|ratio|Ripserer heap|Cubical Ripser heap|
|:-------------|:------|:--|:---------|:-------------|:----|:------------|:------------------|
|`lena512`     |262144 |1  |787 ms    |299 ms        |2.631|145.0 MiB    |49.3 MiB           |
|`lena1999x999`|1997001|1  |2.87 s    |2.009 s       |1.429|514.4 MiB    |186.7 MiB          |
|`bonsai64`    |262144 |2  |2.875 s   |2.996 s       |0.96 |280.6 MiB    |1.3 GiB            |
|`bonsai128`   |2097152|2  |31.151 s  |14.733 s      |2.114|1.5 GiB      |1.9 GiB            |
|`head128`     |2097152|2  |24.102 s  |12.434 s      |1.938|1.5 GiB      |1.9 GiB            |

## Comparison with Eirene

In these benchmarks, we compare Ripserer to
[Eirene.jl](https://github.com/Eetion/Eirene.jl). Ripserer benchmarks were run with
`alg=:involuted`, so this measures the time it takes to compute representative cycles.

|dataset     |size|dim|threshold|Ripserer  |Eirene  |ratio|
|:-----------|:---|:--|:--------|:---------|:-------|:----|
|`gcycle`    |100 |3  |         |6.231 s   |24.158 s|0.258|
|`hiv`       |1088|1  |         |1.824 s   |7.774 s |0.235|
|`dragon1000`|1000|1  |         |575 ms    |8.441 s |0.068|
|`celegans`  |297 |2  |         |4.217 s   |4.588 s |0.919|
|`o3_1024`   |1024|3  |1.8      |5.735 s   |8.314 s |0.69 |
|`random16`  |50  |7  |         |8.577 s   |7.688 s |1.116|
