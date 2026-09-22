#!/bin/sh

open https://github.com/ChristianGaser/CAT-Surface/actions/workflows/binaries.yml
echo "Download last artifacts to ~/Downloads/"
read

[ -f ~/Downloads/Progs ] && rm -r ~/Downloads/Progs

open ~/Downloads/cat-surface-ubuntu-x86_64.zip
sleep 5
rm ~/Downloads/Progs/*.[co]
for i in CAT.glnx86/CAT*; do cp ~/Downloads/Progs/`basename ${i}` CAT.glnx86/ ; done
rm -r ~/Downloads/Progs

open ~/Downloads/cat-surface-macos-x86.zip
sleep 5
rm ~/Downloads/Progs/*.[co]
for i in CAT.maci64/CAT*; do cp ~/Downloads/Progs/`basename ${i}` CAT.maci64/ ; done
rm -r ~/Downloads/Progs

open ~/Downloads/cat-surface-macos-arm64.zip
sleep 5
rm ~/Downloads/Progs/*.[co]
for i in CAT.maca64/CAT*; do cp ~/Downloads/Progs/`basename ${i}` CAT.maca64/ ; done
rm -r ~/Downloads/Progs

open ~/Downloads/cat-surface-windows-x86.zip
sleep 5
rm ~/Downloads/Progs/*.[co]
for i in CAT.w32/CAT*.exe; do cp ~/Downloads/Progs/`basename ${i}` CAT.w32/ ; done
rm -r ~/Downloads/Progs

chmod a+x CAT.*/*

