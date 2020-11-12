#!/bin/bash
# This script generates new gray matter surfaces at indicated fractions of the cortical thickness. Takes the average of white matter surfaces generated from CIVET and an additional surface (e.g. midsurface (generated by CIVET), 25%, or 12.5% surfaces (generated by a previous iteration of this script)

base_surface=$1
old_fraction=$2 #fraction used alongside the white matter to calculate new surface (expressed w/o percentage sign). E.g 25, mid
new_fraction=$3  #new percentile surface fraction of cortical thickness (expressed w/o percentage sign). E.g. 12_5 (corresponds to 25 as old_fraction), 25 (corresponds to mid as old_fraction)

module unload minc-toolkit
module load CIVET/1.1.12

if base_surface==0
then
        cd ./surfaces
        for file in WM_0_surfaces/*left*
        do
                average_surfaces GM_"$new_fraction"_surfaces/$(basename $file _WM_0_surface_left.obj)_GM_"$new_fraction"_surface_left.obj none none 1 $file GM_"$old_fraction"_surfaces/$(basename $file _WM_0_surface_left.obj)_GM_"$old_fraction"_surface_left.obj	
        done
        for file in WM_0_surfaces/*right*
        do
                average_surfaces GM_"$new_fraction"_surfaces/$(basename $file _WM_0_surface_right.obj)_GM_"$new_fraction"_surface_right.obj none none 1 $file GM_"$old_fraction"_surfaces/$(basename $file _WM_0_surface_right.obj)_GM_"$old_fraction"_surface_right.obj
        done

elif base_surface==12_5
then
        cd ./surfaces
        for file in GM_12_5_surfaces/*left*
        do
                average_surfaces GM_"$new_fraction"_surfaces/$(basename $file _GM_12_5_surface_left.obj)_GM_"$new_fraction"_surface_left.obj none none 1 $file GM_"$old_fraction"_surfaces/$(basename $file _WM_0_surface_left.obj)_GM_"$old_fraction"_surface_left.obj	
        done
        for file in GM_12_5_surfaces/*right*
        do
                average_surfaces GM_"$new_fraction"_surfaces/$(basename $file _GM_12_5_surface_right.obj)_GM_"$new_fraction"_surface_right.obj none none 1 $file GM_"$old_fraction"_surfaces/$(basename $file _WM_0_surface_right.obj)_GM_"$old_fraction"_surface_right.obj
        done
fi