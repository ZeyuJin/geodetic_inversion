#!/bin/bash

if [[ $# -ne 3 ]]; then
echo ""
echo "  Usage: $0 SAT master.PRM dem.grd"
echo "  Get satellite look angle components as grd files, based on dem"
echo "  SAT can be ENVI, ALOS, or ERS"
echo "  master.PRM must be found in the SLC directory"
echo "  Output: lle.grd, lln.grd, llu.grd (same dimensions as dem.grd)"
echo ""
exit 1
fi
SAT=$1
PRM=$(basename $2)
dem=$3
##range='-R242.6500/245.7000/32.599997877/34.299999237'
range=`gmt grdinfo -I- $dem`
int=`gmt grdinfo -I $dem`

echo "grd2xyz $dem | eval ${SAT}_look $PRM > look.xyz"
gmt grd2xyz $dem | eval ${SAT}_look $PRM > look.xyz
#awk '{print $1+360,$2,$4}' look.xyz > look.lle
#awk '{print $1+360,$2,$5}' look.xyz > look.lln
#awk '{print $1+360,$2,$6}' look.xyz > look.llu
awk '{print $1,$2,$4}' look.xyz > look.lle
awk '{print $1,$2,$5}' look.xyz > look.lln
awk '{print $1,$2,$6}' look.xyz > look.llu
echo "blockmedian look.lle $range $int -bo | surface -bi -Glke.grd $int $range -T0.5 -V"
gmt blockmedian look.lle $range $int -bo | gmt surface -bi -Glke.grd $int $range -T0.5 -V
gmt blockmedian look.lln $range $int -bo | gmt surface -bi -Glkn.grd $int $range -T0.5 -V
gmt blockmedian look.llu $range $int -bo | gmt surface -bi -Glku.grd $int $range -T0.5 -V

rm -r look.*

#add unit look vectors to the set of observation points:
rm -f gps_los.dat
gmt grdtrack -Glke.grd -Glkn.grd -Glku.grd gps_sel.dat > gps_los.dat

# format for Matlab (single space delimeters):
circ gps_los.dat junk 1 7
mv junk gps_los.dat

echo ""
echo " Finished!"
