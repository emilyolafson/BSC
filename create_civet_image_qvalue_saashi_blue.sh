#!/bin/bash
#USAGE
#create_civet_image.sh left.obj left_statmap left_column left_FDR right.obj right_statmap right_column right_FDR output.png

#Automatically generates a standardized view of a brain, using stats files for left and right hemispheres (from RMINC)
#Colourizes objects using stats files, ray traces standard views and them merges them together

set -euo pipefail
#IFS=$'\n\t'

err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

awkscript='BEGIN{
    OFS=FS=","
    split(target, t_targets, ",")
    for (i in t_targets)
        targets[t_targets[i]] = i
}
NR==1 {
    for (i = 1; i <= NF; i++) head[i] = $i
}

NR !=1{
    for (i = 1; i <= NF; i++) {
        if (head[i] in targets){
            print $i
        }
    }
}'

if [[ $# -ne 9 ]]
then
echo Usage: create_civet_image_qvalue_saashi.sh left_object.obj left_statmap left_column left_lower_thresh,left_upper_thresh right_object.obj right_statmap right_column right_lower_thresh,right_upper_thresh PREFIX
exit
fi


TMPDIR=$(mktemp -d)
leftbrain=$1
leftstatmap=$2
leftcolumn=$3
leftFDR5=$(echo $4 | cut -d"," -f 1)
leftFDR1=$(echo $4 | cut -d"," -f 2)
rightbrain=$5
rightstatmap=$6
rightcolumn=$7
rightFDR5=$(echo $8 | cut -d"," -f 1)
rightFDR1=$(echo $8 | cut -d"," -f 2)
prefix=$9

#Check file type
filename=$(basename $leftstatmap)
extension="${filename##*.}"

#Pick out column, multiply by mask and store in new tempfile
case $extension in
    csv)
        cut -d "," -f $leftcolumn $leftstatmap | tail -n +2 | paste -d\+ /opt/quarantine/resources/CIVET/CIVET-CC-mask-inv.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $leftstatmap)
        cut -d "," -f $rightcolumn $rightstatmap | tail -n +2 | paste -d\+ /opt/quarantine/resources/CIVET/CIVET-CC-mask-inv.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $rightstatmap)
        ;;
    vertstats)
        cut -d " " -f $leftcolumn $leftstatmap | tail -n +4 | paste -d\+ /opt/quarantine/resources/CIVET/CIVET-CC-mask-inv.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $leftstatmap)
        cut -d " " -f $rightcolumn $rightstatmap | tail -n +4 | paste -d\+ /opt/quarantine/resources/CIVET/CIVET-CC-mask-inv.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $rightstatmap)
        ;;
    txt)
        cut -d " " -f $leftcolumn $leftstatmap |  paste -d\+ /opt/quarantine/resources/CIVET/CIVET-CC-mask-inv.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $leftstatmap)
        cut -d " " -f $rightcolumn $rightstatmap | paste -d\+ /opt/quarantine/resources/CIVET/CIVET-CC-mask-inv.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $rightstatmap)
        ;;
esac


#Redefine statmap to masked tempfile
leftstatmap=$TMPDIR/$(basename $leftstatmap)
rightstatmap=$TMPDIR/$(basename $rightstatmap)


#Colourize left
colour_object $leftbrain $leftstatmap $TMPDIR/left_comb.obj user /home/cic/olaemi/emilyO/Tissue_contrast/analysis/main_effect/smoothed_10mm/virdis.lut $leftFDR1 $leftFDR5 "0 1 1" white replace


#Colourize right
colour_object $rightbrain $rightstatmap $TMPDIR/right_comb.obj red $rightFDR1 $rightFDR5 red white replace


#Generate views
echo """ray_trace -output $TMPDIR/left_medial.rgb $TMPDIR/left_comb.obj -size 1000 1000 -crop -sup 3 -bg '1 0 1' -right -shadows -directional 1 -1 -1 1 1 1 -directional -1 -1 -1 1 1 1 -directional -1 -1 1 1 1 1
ray_trace -output $TMPDIR/left_lateral.rgb $TMPDIR/left_comb.obj -size 1000 1000 -crop -sup 3 -bg '1 0 1' -left -shadows -directional 1 -1 -1 1 1 1 -directional -1 -1 -1 1 1 1 -directional -1 -1 1 1 1 1
ray_trace -output $TMPDIR/right_medial.rgb $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg '1 0 1' -left -shadows -directional 1 -1 -1 1 1 1 -directional -1 -1 -1 1 1 1 -directional -1 -1 1 1 1 1
ray_trace -output $TMPDIR/right_lateral.rgb $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg '1 0 1'  -right -shadows -directional 1 -1 -1 1 1 1 -directional -1 -1 -1 1 1 1 -directional -1 -1 1 1 1 1""" | parallel

#Create pngs
for file in $TMPDIR/*.rgb
do
    echo convert $file -fill transparent -fuzz 25\% -draw \'color 10,10 floodfill\' $TMPDIR/$(basename $file rgb)png
done | parallel

for file in $TMPDIR/*png
do
    cp -f $file ${prefix}_$(basename $file)
done

rm -rf $TMPDIR
