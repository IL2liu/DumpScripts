#/bin/sh

################################################################
# map the NA from gem mapping into binary profile by HMM3state #
################################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

for f in $(ls *.bin)
do
   # get file without extension, i.e, 109 from 109.txt
   f_base=${f%%.*}
   map_f=$f_base".na"
   # combine binary profile with profile that has NA
   paste -d'\t' $f $map_f  | awk -F'\t' 'BEGIN { OFS="\t" } {print $1,$2,$3,$4,$8}' > tmp
   awk 'BEGIN{ OFS="\t" }{if ($5 == "NA") print $1,$2,$3,$4="NA"; else print $1,$2,$3,$4}' tmp > tmp2
   map_f=$f_base".bed"
   mv tmp2 $map_f
done
