#!/bin/bash 

##
# @Description: 
# Rearrange result files from jatos
# such that they are listed under subject id
#   sourcedir - highest level dir of the study results
#   destdir - (noexisting) directory for the files
##

# default values
DESTDIR="./"

function show_usage() {
    printf "Usage: $0 sourcedir destdir\n"
    #printf "\n"
    #printf "Options:\n"
    #printf " -s|--sourcedir, source directory\n"
    #printf " -d|--destdir, destination directory\n"
    #printf " -h|--help, Print help\n"

    return 0
}

if [ -z "$1" ]; then # if parameter 1 is empty
    show_usage
    exit 0
else
    SOURCEDIR=$1
    if [ -n "$2" ]; then # if parameter 2 is not empty
        DESTDIR=$2
    fi
    ids=($(find $SOURCEDIR -type f | egrep -o '/[a-zA-Z0-9]+_' | egrep -o '[a-zA-Z0-9]+' | sort | uniq))
fi

for id in "${ids[@]}"
do
    echo "Subject: $id"
    days=($(find $SOURCEDIR -type f | egrep -o "${id}_day0[0-3]" | uniq))
    ndays=${#days[@]}

    for day in `seq 1 $ndays`
    do
        dirname="${DESTDIR}/${id}/day0${day}"
        if [ ! -d $dirname ]; then
            mkdir -p $dirname
            
            files=($(find $SOURCEDIR -regex ".*/${id}_day0${day}.*"))
            for file in "${files[@]}"
            do
                filename="${file##*/}"
                out="${id}_day0${day}_"
                name="${filename//${out}/}"
                nfile="${dirname}/${name}"
                # check if the file already exists
                if [ -f "$nfile" ]; then
                    e=($(date +%s))
                    echo "A file ${name} already exists"
                    cp $file "${dirname}/_${e}.${name}"
                else 
                    cp $file "${dirname}/${name}"
                fi 
                
            done

            # copy file with pictures to nodes assignment
            if [ "$day" == "1" ]; then
                
                file="$(find $SOURCEDIR -regex ".*/${id}_pics2nodesObj.json")"
                filename="${file##*/}"
                out="${id}_"
                name="${filename//${out}/}"

                cp $file "${dirname}/${name}"
            fi      
        fi 
    done
done



