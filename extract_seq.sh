#!/bin/bash
out=$1 # output directory that contains the missing barcode files
barcode_file=$2
inp_read1=$3
inp_read2=$4

f1=${out}/missing_barcode_piscem.txt
f2=${out}/missing_barcode_chromap.txt

read1_piscem=${out}/read1_miss_piscem.fq
barcode_piscem=${out}/barcode_miss_piscem.fq
read2_piscem=${out}/read2_miss_piscem.fq ## will be used only if inp_read2 is provided

read1_chromap=${out}/read1_miss_chromap.fq
barcode_chromap=${out}/barcode_miss_chromap.fq
read2_chromap=${out}/read2_miss_chromap.fq ## will be used only if inp_read2 is provided

headers_missing_chromap=${out}/headers_missing_chromap.txt

printf "" > ${read1_piscem}
printf "" > ${barcode_piscem}
if [[ -n $4 ]]
    then
        printf "" > ${read2_piscem}
fi

### extract each barcode and find its corresponding header in the barcode file
while IFS= read -r query; do
    # Use grep to search for the query in data.txt
    # header_barcode=$(grep -B 1 $query ${barcode_file} | grep -v -- -- | awk 'NR %2 != 0')
    
    # header_f1=$(echo ${header_barcode} | sed 's/ 2:N:0/ 1:N:0/')
    # grep -h -A 3 "${header_f1}" ${inp_read1} >> ${read1_piscem}
    # grep -h -A 3 "${header_barcode}" ${barcode_file} >> ${barcode_piscem}
    # if [[ -n $4 ]]
    #     then
    #         header_f2=$(echo ${header_barcode} | sed 's/ 2:N:0/ 3:N:0/')
    #         grep -h -A 3 "${header_f2}" ${inp_read2} >> ${read2_piscem}

    # fi
    mapfile -t header_barcodes < <(grep -B 1 $query ${barcode_file} | grep -v -- -- | awk 'NR %2 != 0')
    count_barcode=${#header_barcodes[@]}
    if [ $count_barcode -gt 1 ]; then
        echo $count_barcode
        echo "barcode $query repeated ${count_barcode} times" >> ${out}/barcode_repeat_piscem.txt
    fi
    for header_barcode in "${header_barcodes[@]}";do
        header_f1=$(echo ${header_barcode} | sed 's/ 2:N:0/ 1:N:0/')
        grep -h -A 3 "${header_f1}" ${inp_read1} >> ${read1_piscem}
        grep -h -A 3 "${header_barcode}" ${barcode_file} >> ${barcode_piscem}
        if [[ -n $4 ]]
            then
                header_f2=$(echo ${header_barcode} | sed 's/ 2:N:0/ 3:N:0/')
                grep -h -A 3 "${header_f2}" ${inp_read2} >> ${read2_piscem}

        fi
    done
done < $f1

printf "" > ${read1_chromap}
printf "" > ${barcode_chromap}
if [[ -n $4 ]]
    then
        printf "" > ${read2_chromap}
fi

while IFS= read -r query; do
    # Use grep to search for the query in data.txt
    mapfile -t header_barcodes < <(grep -B 1 $query ${barcode_file} | grep -v -- -- | awk 'NR %2 != 0')
    count_barcode=${#header_barcodes[@]}
    if [ $count_barcode -gt 1 ]; then
        echo $count_barcode
        echo "barcode $query repeated ${count_barcode} times" >> ${out}/barcode_repeat_chromap.txt
    fi
    for header_barcode in "${header_barcodes[@]}";do
        header_f1=$(echo ${header_barcode} | sed 's/ 2:N:0/ 1:N:0/')
        grep -h -A 3 "${header_f1}" ${inp_read1} >> ${read1_chromap}
        grep -h -A 3 "${header_barcode}" ${barcode_file} >> ${barcode_chromap}
        if [[ -n $4 ]]
            then
                header_f2=$(echo ${header_barcode} | sed 's/ 2:N:0/ 3:N:0/')
                grep -h -A 3 "${header_f2}" ${inp_read2} >> ${read2_chromap}

        fi
    done
done < $f2
