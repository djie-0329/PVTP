# Parallel_and_Approximate_VTP_source_code

The source code is a reference implementation of the following paper:
Jie Du, Ying He, Zheng Fang, Wenlong Meng, Shiqing Xin. On Vertex-oriented Triangle Propagation Algorithm: Parallelization and Approximation. Under review. 2020.

## Usage:
-alg [algorithmIndex]: the index of algorithm, 0: VTP, 1: Parallel-VTP (PVTP), 2: Approximate-VTP (AVTP), 3: Parallel Approximate-VTP (PAVTP).

-m [meshFIle]: input model file (only support .obj files).

-s [sourceIndex]: the index of source.

-np [numProcs]: number of threads.

-l [lambda]: parameter for Approximate VTP algorithm.

-o [output]: bool: to output geodesic distance result.

## Examples:
```
VTP:          -alg 0 -m bunny_nf144k.obj -s 0 -o
```
```
PVTP:         -alg 1 -m bunny_nf144k.obj -s 0 -np 4 -o
```
```
AVTP:         -alg 2 -m bunny_nf144k.obj -s 0 -l 8 -o
```
```
PAVTP:        -alg 3 -m bunny_nf144k.obj -s 0 -l 8 -np 4 -o
```

## Code Instructions

```
geodesic_algorithm_exact.h: original VTP
```
```
geodesic_algorithm_parallel_fwp_exact.h: parallel VTP (PVTP)
```
```
geodesic_algorithm_approximate.h: approximate VTP (AVTP)
```
```
geodesic_algorithm_parallel_fwp_approximate.h: parallel and approximate VTP (PAVTP)
```

## More

Our parallel implementation utilizes the tbb library "tbb2019_20190605oss".

Please cite this paper if you use this code:
```
Jie Du, Ying He, Zheng Fang, Wenlong Meng, Shiqing Xin. On Vertex-oriented Triangle Propagation Algorithm: Parallelization and Approximation. Under review. 2020.
```


