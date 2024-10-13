Efficient Cones
===

This is a C++ implementation of the cone singularity construction in the following paper:

[Efficient Cone Singularity Construction for Conformal Parameterizations](https://doi.org/10.1145/3618407).
Mo Li, [Qing Fang](https://qingfang1208.github.io/), Zheng Zhang, [Ligang Liu](http://staff.ustc.edu.cn/~lgliu/), [Xiao-Ming Fu](https://ustc-gcl-f.github.io/).
*ACM Transactions on Graphics (SIGGRAPH Asia)*, 42(6), 2023.

## Dependencies
* [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/)
* [Eigen](http://eigen.tuxfamily.org/)
* [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)

## Quick Start

```bash
git clone --recursive git@github.com:pielet/MovingCones.git
```
After cloning the repo, build and install cholmod seperatly. Then

```bash
mkdir build && cd build
cmake ..
make
./MovingCones [INPUT_OBJ] [OUTPUT_PATH] [INPUT_DISTORTION_BOUND]
```
Command arguments:
* [INPUT_OBJ]: a input mesh (.obj).
* [OUTPUT_PATH]: the output path.
* [INPUT_DISTORTION_BOUND]: the input area distortion bound (default: 0.2).

## Output
* cones.txt: the first column records cones ids(index start from 1) and the second column records the non-zero Gaussian curvature.