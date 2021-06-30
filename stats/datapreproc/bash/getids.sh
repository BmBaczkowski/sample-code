#!/bin/bash 

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

filelistname="$SOURCEDIR/list.txt"
arr=()
if [ -f $filelistname ]; then
while IFS= read -r line; do
   arr+=("$line")
done < $filelistname
fi

# find new id that is not yet in the list 
narr=(`echo ${ids[@]} ${arr[@]} | tr ' ' '\n' | sort | uniq -u`)
if [ ${#narr[@]} -eq 0 ]; then
    echo "There are no new ids. The life just goes on."
else
    echo "There are some new ids. I will append the list."
    for id in "${narr[@]}"
    do echo $id >> $filelistname
    done
fi