# Benchmarks

The following tables show benchmarks that compare Ripserer's performance with
[Ripser](https://github.com/Ripser/ripser) and [Cubical
Ripser](https://github.com/CubicalRipser/). The benchmarking code and more info about the
datasets is available [here](https://github.com/mtsch/RipsererBenchmarks.jl).

All benchmarks were performed on a laptop with an Intel(R) Core(TM) i5-4200U CPU @ 1.60GHz
with 8GB of RAM.

We used [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl/) to perform the
timing benchmarks and [Valgrind's Massif
tool](https://www.valgrind.org/docs/manual/ms-manual.html) to measure peak heap sizes
(i.e. total memory footprint).

All benchmarks include I/O, so there might be some noise in the smaller benchmarks. For
Ripserer, [CSV.jl](https://github.com/JuliaData/CSV.jl) was used to read the
files. Ripserer's memory footprint benchmarks include the Julia runtime, which tends to use
around 100MiB memory.

The benchmarks were performed with Ripserer v0.14.3 and `master` versions of Ripser (commit
hash `286d369`) and Cubical Ripser (commit hashes `6edb9c5` for 2D and `a063dac` for 3D).

## Vietoris-Rips Persistent Homology

In this experiment, we performed benchmarks with the datasets presented in the [Ripser
article](https://arxiv.org/abs/1908.02518). We only used the datasets that we were able to
run with less than 8GB memory. `ripserer_s` denotes `ripserer` run with `sparse=true`. All
datasets were parsed as `Float32` as that is what Ripser supports.

The timing results are in the table below. The ratio column shows the better ratio.

| dataset    | size  | `dim_max` | `threshold` | `ripserer` | `ripserer_s` | `ripser`  | ratio |
| :--------- | :---- | :-------- | :---------- | :--------- | :----------- | :-------- | :---- |
| `random16` | 50    | 2         |             | 8.1ms      | 9.3ms        | 10.2 ms   | 0.80  |
| `o3_1024`  | 1024  | 3         | 1.8         | 4.3s       | 3.1s         | 3.1 s     | 1.01  |
| `o3_4096`  | 4096  | 3         | 1.4         | 143.2s     | 76.3s        | 75.6 s    | 1.01  |
| `fract-r`  | 512   | 2         |             | 22.4s      | 25.4s        | 19.6 s    | 1.14  |
| `sphere_3` | 192   | 2         |             | 1.6s       | 1.8s         | 1.5 s     | 1.04  |
| `dragon`   | 2000  | 1         |             | 3.2s       | 3.8s         | 2.8 s     | 1.12  |

The next table shows peak heap sizes as reported by `valgrind`.

| dataset    | size | `dim_max` | `threshold` | `ripserer` | `ripserer_s` | `ripser`  |
| :--------- | :--- | :-------- | :---------- | :--------- | :----------- | :-------- |
| `random16` | 50   | 2         |             | 111.1 MiB  | 111.1 MiB    | 1.1 MiB   |
| `o3_1024`  | 1024 | 3         | 1.8         | 374.1 MiB  | 418.2 MiB    | 151.0 MiB |
| `o3_4096`  | 4096 | 3         | 1.4         | 4.7 GiB    | 4.9 GiB      | 4.1 GiB   |
| `fract-r`  | 512  | 2         |             | 2.2 GiB    | 2.2 GiB      | 2.0 GiB   |
| `sphere_3` | 192  | 2         |             | 287.0 MiB  | 288.5 MiB    | 209.5 MiB |
| `dragon`   | 2000 | 1         |             | 316.7 MiB  | 350.4 MiB    | 296.8 MiB |

## Alpha-Rips Persistent Homology

Here, we used the 1-skeleton of the Delaunay triangulation of a point cloud as a sparse
distance matrix. This benchmark is intended to benchmark performance in very sparse cases,
as triangulations only have ``\mathcal{O}(n)`` ``d``-simplices. All datasets were parsed as
`Float32` as that is what Ripser supports.

Timing results:

| dataset    | size  | `dim_max` | `ripserer` | `ripser` | ratio |
| :--------- | :---- | :--------- | :-------- | :------- | :---- |
| `torus`    | 10000 | 2         | 937ms      | 1.2s     | 0.79  |
| `3_sphere` | 3000  | 3         | 623ms      | 809ms    | 0.77  |
| `4_sphere` | 2000  | 4         | 5.9s       | 6.3s     | 0.94  |
| `5_sphere` | 1000  | 5         | 49.7s      | 47.6s    | 1.04  |
| `dragon`   | 2000  | 2         | 58ms       | 76.3s    | 0.76  |

Peak heap size:

| dataset    | size  | `dim_max` | `ripserer` | `ripser`  |
| :--------- | :---- | :-------- | :--------- | :-------- |
| `torus`    | 10000 | 2         | 138.4 MiB  | 33.2 MiB  |
| `3_sphere` | 3000  | 3         | 130.0 MiB  | 27.7 MiB  |
| `4_sphere` | 2000  | 4         | 387.2 MiB  | 202.0 MiB |
| `5_sphere` | 1000  | 5         | 2.4 GiB    | 1.5 GiB   |
| `dragon`   | 2000  | 2         | 110.9 MiB  | 33.2 MiB  |

## Cubical Persistent Homology

In this benchmark, we use some of the datasets presented in the [Cubical
Ripser](https://arxiv.org/abs/2005.12692) article. We limited the 2D image size to 1999×999
as the current `master` (commit hash `6edb9c5`) version of 2D Cubical Ripser throws an
assertion error for anything larger. We were also unable to perform 3D 256×256×256 image
benchmarks due to Ripserer running out of memory. The `eltype` of all datasets is `Float64`,
because that is what Cubical Ripser supports. When running Ripserer in the real world, it's
a good idea to use the native data types. This will slightly reduce the memory footprint and
increase performance.

Timing results:

| dataset  | size        | `ripserer` | `cubical_ripser` | ratio |
| :------- | :---------- | :-------   | :--------------- | :---- |
| `lena`   | 1999×999    | 3.8 s      | 2.0 s            | 1.90  |
| `bonsai` | 128×128×128 | 33.1 s     | 15.0 s           | 2.21  |
| `lena`   | 512×512     | 842.6 ms   | 295.9 ms         | 2.85  |
| `head`   | 128×128×128 | 26.0 s     | 12.8 s           | 2.03  |
| `bonsai` | 64×64×64    | 3.0 s      | 3.0 s            | 1.01  |

Peak heap size:

| dataset  | size        | `ripserer` | `cubical_ripser` |
| :------- | :---------- | :--------- | :--------------- |
| `lena`   | 512×512     | 145.0 MiB  | 49.3 MiB         |
| `lena`   | 1999×999    | 514.4 MiB  | 186.7 MiB        |
| `bonsai` | 64×64×64    | 280.6 MiB  | 1.3 GiB          |
| `bonsai` | 128×128×128 | 1.5 GiB    | 1.9 GiB          |
| `head`   | 128×128×128 | 1.5 GiB    | 1.9 GiB          |
