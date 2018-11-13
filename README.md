# CoorTransform
Used to change wgs84/cgcs2000 to gauss coordinates

# Depends
It depends on integrate lib(gsl).

# Usage
See the \.cc file

# Run
clang++ coor_transform.cc -std=c++17 -I <GSL header> -L <GSL Lib> -lgsl -lgslcblas
