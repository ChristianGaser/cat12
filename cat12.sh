#! /bin/bash
# ______________________________________________________________________
#
# Christian Gaser, Robert Dahnke
# Structural Brain Mapping Group (https://neuro-jena.github.io)
# Departments of Neurology and Psychiatry
# Jena University Hospital
# ______________________________________________________________________
# $Id$
version='cat12.sh $Id$'

echo "##############################################################"
echo "   cat12.sh is deprecated. Please now use cat_batch_cat.sh.   "
echo "##############################################################"

cat12_dir=$(dirname "$0")
args=("$@")

${cat12_dir}/cat_batch_cat.sh ${args}

