#!/bin/bash

## good single launch:
./locrun 1 144 1976 0.0022 > loc.dat

gnuplot distance.gnu
gnuplot loc.gnu
	
