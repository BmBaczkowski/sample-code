#! /opt/miniconda3/bin/python

import sys, os, fnmatch, re
import json
from shutil import copyfile
from datetime import datetime

def show_usage(name):
    print ('Usage: ', name, 'sourcedir destdir\n')
    return sys.exit()

# default values
DESTDIR="./results"

# read out parameters
n = len(sys.argv) 
if (n == 1):
    name = sys.argv[0]
    show_usage(os.path.basename(name))
else:
    SOURCEDIR = sys.argv[1]
    if (n > 2):
        DESTDIR = sys.argv[2]
    if not os.path.isdir(DESTDIR):
        print('Creating destdir...')
        os.mkdir(DESTDIR)
    jsonkeyfilename = os.path.join(DESTDIR, 'key.json')

# find all ids from the available JATOS files
print("Generating a list of subject ids...")
ids = set()
for root, dirs, files in os.walk(SOURCEDIR):
    for f in fnmatch.filter(files, '*_pics2nodes.json'): # should make it faster since there is only 1 file
        match = re.findall("[a-zA-Z0-9]+_", f)
        pid = re.findall("[a-zA-Z0-9]+", match[0])[0]
        ids.add(pid)
# convert ids to a list
ids = list(ids)

# read file of existing ids
if os.path.isfile(jsonkeyfilename):
    with open(jsonkeyfilename, 'r') as myfile:
        data = myfile.read()
        mydic = json.loads(data) 
    keys = list(mydic.keys())
    newids = set(ids)-set(keys)
    newids = list(newids)
    if bool(newids):
        vals = list(mydic.values())
        vals = [int(x) for x in vals]
        n = [str(i).zfill(3) for i in range(max(vals)+1,max(vals)+1+len(newids)+1)]
        nmydic = dict(zip(newids, n))
        mydic.update(nmydic) 
    
else:
    n = [str(i).zfill(3) for i in range(1,len(ids)+1)]
    mydic = dict(zip(ids, n))


# save id key as a JSON object
print("Saving the list of ids in a key.json file...")
with open(jsonkeyfilename, 'w') as json_file:
    json.dump(mydic, json_file)

# move files
print("Moving files...")
for subject_id_orig, subject_id_new in mydic.items():
    print("Subject: ", subject_id_orig, "->", subject_id_new)

    dirname="%s/%s/" % (DESTDIR, subject_id_new)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
        basename = "%s_" % (subject_id_orig)
        # find all files
        for root, dirs, files in os.walk(SOURCEDIR):
            for f in fnmatch.filter(files, "%s_*" % (subject_id_orig)):
                x = re.sub(basename, '', f)
                sourcefile = os.path.join(root, f)
                destfile = os.path.join(dirname, x)
                if not os.path.isfile(destfile):
                    copyfile(sourcefile, destfile)
                else:
                    # filename based on file modification date
                    destfile = os.path.join(dirname, "%i_%s" % (os.path.getmtime(sourcefile), x))
                    print("File already exists")
                    copyfile(sourcefile, destfile)
    else:
        print("Data from this subject and session already exist")
