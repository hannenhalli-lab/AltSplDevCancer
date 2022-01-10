#!/bin/bash


cancerType=$1
indexFile=$2
ProjectCode=$3
allocatedCPUs=$4
URL=$5

pathTodependencies="/home/singha30/developmentalDataProcess/kallistoRequirements"
mkdir -p /data/singha30/developmentalDataProcess/RNAseqAnalysis/${cancerType}
pathToCancerType=/data/singha30/developmentalDataProcess/RNAseqAnalysis/${cancerType}

CPUS_TrimGallore=$(($allocatedCPUs / 2))
CPUS_Kallisto=$(($allocatedCPUs))

#echo $ProjectCode
#gdc-client download $ProjectCode -t $pathTodependencies/$tokenFile --config $pathTodependencies/gdcConfigurationFilefor${cancerType}.dtt

pathToDownloadedFolder=${pathToCancerType}/${ProjectCode}
wget -P $pathToDownloadedFolder $URL
TarFile=$(ls $pathToDownloadedFolder | grep "fastq.gz")
gunzip $pathToDownloadedFolder/$TarFile

return=$? ###### if $return is 0, then previous command was successfull and it would be safe to delete pathToTar folder below#####
if [ $return==0 ]
then
	rm $pathToDownloadedFolder/$TarFile
else
	echo "some error happened with decompression output for ${ProjectCode} : program is exiting..."
	break
fi

FASTA=$(ls $pathToDownloadedFolder | grep ".*.fastq$")

trim_galore $pathToDownloadedFolder/$FASTA --cores $CPUS_TrimGallore -o  $pathToDownloadedFolder

return=$?
if [ $return==0 ]
then
	rm $pathToDownloadedFolder/$FASTA
	rm $pathToDownloadedFolder/*trimming_report.txt
else
	echo "some error happened with TrimGalore output for ${ProjectCode} : program is exiting..."
	break
fi
#rm pathToTar #####this will remove all the downloaded data folder content

FASTA=$(ls $pathToDownloadedFolder | grep ".*.fq$")

echo "Now doing mapping to the reference transcriptome with Salmon for ${ProjectCode}"
kallisto quant -i $pathTodependencies/$indexFile -o $pathToDownloadedFolder/Measurements $pathToDownloadedFolder/$FASTA --rf-stranded --bias --single -l 200 -s 30  -t $CPUS_Kallisto 

if [ $return==0 ]
then
	rm $pathToDownloadedFolder/$FASTA
else
	echo "some error happened with salmon output for ${ProjectCode} : program is exiting..."
	break
fi
tar -czvf $pathToCancerType/MeasurementsFor$ProjectCode.tar.gz $pathToDownloadedFolder/Measurements/
#gzip -cvf $pathToCancerType/Measurements/quant.sf > $pathToCancerType/MeasurementsFor_$ProjectCode.tar.gz 
return=$?
if [ $return==0 ]
then
	rm -r $pathToDownloadedFolder
else
	echo "some error happened with compression for ${ProjectCode} : program is exiting..."
	break
fi


