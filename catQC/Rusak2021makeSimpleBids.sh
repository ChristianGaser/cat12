# prepare Rusak

maindir=20211122-SyntheticDataset
resdir=Rusak2021/SimAthrophy
subjects=$(find $maindir  -name "*_S_*" -depth 1 )
 
for S in $subjects
do
	SN=$(basename $S | tr -d _)
	SN2=sub-ADNI$SN

	mkdir -p $resdir/$SN2

	SES=$(ls $S) 
	for SI in $SES
	do
		SN3=$(printf ${SN2}_simGMatrophy%0.2fmm.nii.gz $SI )
		cp  ./$S/$SI/SyntheticMRI.nii.gz  $resdir/$SN2/$SN3
	done
	
	gunzip $resdir/$SN2/*.nii.gz
	exit
done

