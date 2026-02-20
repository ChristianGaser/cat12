#! /bin/bash
# Tool for adding revision and date to m-files that were prepared in CAT.
# ______________________________________________________________________
#
# Christian Gaser, Robert Dahnke
# Structural Brain Mapping Group (https://neuro-jena.github.io)
# Departments of Neurology and Psychiatry
# Jena University Hospital
# ______________________________________________________________________


REVISION=`git rev-list --count HEAD`
DATE=`git log --date short |grep "Date:"|head -1|cut -f2 -d':'|sed -e s'/ //g'`

perl -p -i -e 's/\$Id\$/\$Id: '$REVISION' '$DATE' \$/g' CAT/*.m
perl -p -i -e 's/\$Id\$/\$Id: '$REVISION' '$DATE' \$/g' CAT/*.sh