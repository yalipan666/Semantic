#!/bin/bash

# set parameters
declare -a tSubjects=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34") # subject index  "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44" "45" "46" "47" "48" "49" "50" "51" 
f="Fixation_Related_Fields_sv"	# your MATLAB function name
t="00:09:99" 	# time: hh-mm-ss
m="60GB" 			# memory

# do not change code below this line
fname="$f.sbatch" 
echo "#!/bin/bash" > $fname
echo "#SBATCH --ntasks 1" >> $fname
echo "#SBATCH --time $t" >> $fname
echo "#SBATCH --qos bbdefault" >> $fname # bbdefault bbshort<10mins
echo "#SBATCH --mem $m" >> $fname
echo "" >> $fname
echo "set -e" >> $fname
echo "" >> $fname
echo "module purge; module load bluebear" >> $fname
echo "module load MATLAB/2020a" >> $fname
echo "" >> $fname
echo "i=\$1" >> $fname
echo "i=\"\${i:1:\${#i}}\"" >> $fname
echo "" >> $fname
echo "matlab -nodisplay -r \"$f(\$i)\"" >> $fname

for i in "${tSubjects[@]}" 
do
	sbatch "$fname" ="$i"
done
