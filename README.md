# Polyhedral-Net Splines for Iso-Geometric Analysis #

C++ library for efficiently solving elliptic partial differential equations on smooth free-form surfaces represented by polyhedral control nets with irregular configurations (grid, n-gon, polar, etc.). For technical details regarding the codebase, please refer to our journal paper [Polyhedral control-net splines for analysis](https://www.sciencedirect.com/science/article/pii/S0898122123004261). This codebase was designed through the contributions of Bhaskar Mishra and William Gregory, as well as the guidance of Kyle Lo and Dr. Jorg Peters.

# Environment

### Recommended Operating Systems
Linux (Ubuntu 20.04 LTS) & macOS (Monterey 12.6.2)
 & Windows 10

### Tools & Dependencies

[CMake](https://cmake.org/) (>=3.9)
[OpenMesh](<https://www.openmesh.org/>) (8.1)
[Eigen](eigen.tuxfamily.org) (3.4)

**Note:**

* CMake will automatically download and install both OpenMesh and Eigen into `/Source/External` using the source with commit hash pointing to the tested version.
* The program should be compatible to equivalent or higher versions with little or no modification.

### Tested Environments
macOS: Apple clang++-11

Check: ubuntu 18.04: g++-6 g++-7 g++-8

Check: ubuntu 20.04: g++-9 g++-10 g++-11

<!--
Check:  Debian 9: g++-6

Check: Debian 10: g++-7 g++-8

Check: Debian 11: g++-9 g++-10 g++11

Check: CentOS 8: g++-8
-->

Check: Windows 10: Visual Studio 2017

# Building #

For UNIX-based system:
```shell
git clone https://[username]@bitbucket.org:surflab/23psiga.git
cd 23psiga
mkdir build
cd build
cmake ../Source
make
```

For Windows:

1. Launch `x86 Native Tools Command Prompt for VS 2017` with admin

2. Run the following commands

```
git clone https://[username]@bitbucket.org:surflab/23psiga.git
cd 23psiga
mkdir build
cd build
cmake ../Source
```

<!--**TODO: Prior to making this public, update above git clone steps to clone this repository rather than the base polyhedral splines**-->

3. In the build folder find and launch `PolyhedralSplinesIGA.sln` with Visual Studio 2017 (or just double click on .sln file)

4. In Visual Studio 2017, set configuration to `Release` mode and switch platform to `Win32` then click `Build Solution`

5. `cd ~/build/` (WARNING: calling executable from other directory might fail.)


**Note:** macOS users need to make sure $PATH includes path to qt5 bin folder

**Note:** alternatively in git: clone HTTPS  cut & paste in terminal;     generate passwd:   gear, personal, (left) app passwords, my label, repository: rw admin

# Usage #
* Input:  quad-dominant mesh in .obj file format
* Output: BB-coefficients written in .bv file format.


Please find [BView file introduction](https://www.cise.ufl.edu/research/SurfLab/bview/#file-format) on UF CISE SurfLab website

## Execution
For UNIX-based system:
```shell
./PolyhedralSplinesIGA /path/to/filename.obj
```
Note: test .obj files are in `/testfile`.

User can add option `-d` or `--DEGREE_RAISE` to raise the degree of all patches to a uniform 3x3
```shell
./PolyhedralSplinesIGA -d /path/to/filename.obj
```

For Windows:
```
Release\PolyhedralSplinesIGA.exe \path\to\filename.obj
```

```
Release\PolyhedralSplinesIGA.exe -d \path\to\filename.obj
```

## View .sf file

copy to sf.dat  in MATLAB:

```shell
load sf.dat;  [r,c] = size(sf);
N = 11; %  10 x 10 quads per patch
N2 = N^2; pats = r/N2  % number of patches
idx = [1:N2]; tmpl = zeros(N,N);
clf;
for jj=1:pats,
   for ii=1:4,
      x{ii}{jj} = tmpl;
      x{ii}{jj}(:) = sf((jj-1)* N2+idx,ii);
   end;
   hsf = surf(x{1}{jj}, x{2}{jj}, x{3}{jj}, x{4}{jj});
   hold on;
end;
axis equal; colormap(hot); shading interp
```


## View .bv file
Users can display .bv files by using either the online or the desktop version of [BView](https://www.cise.ufl.edu/research/SurfLab/bview/).
BView provides surface inspection tools such as zoom, groupwise patch color, Gaussian curvature, highlight lines and more.

Quick view of .bv file:

1. Go to [WebGL](https://www.cise.ufl.edu/research/SurfLab/bview/webgl3_new_2/)

2. Click on `Load File` > `Choose File` to import .bv file

3. Set `Display` > `Patch Detail` to 16x16
