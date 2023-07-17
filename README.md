# CAT: Computational Anatomy Toolbox
This toolbox is a an extension to [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Wellcome Department of Cognitive Neurology) to provide computational anatomy. This covers diverse morphometric methods such as voxel-based morphometry (VBM), surface-based morphometry (SBM), deformation-based morphometry (DBM), and region- or label-based morphometry (RBM).

It is developed by Christian Gaser and Robert Dahnke (Jena University Hospital, Departments of Psychiatry and Neurology) and free but copyright software, distributed under the terms of the [GNU General Public Licence](http://www.gnu.org/licenses/gpl-2.0.html) as published by the Free Software Foundation; either version 2 of the Licence, or (at your option) any later version.

## Download
[CAT12 toolbox](http://141.35.69.218/cat12/cat12_latest.zip)

[CAT12 manual](http://141.35.69.218/cat12/CAT12-Manual.pdf)

[Quick Start Guide](https://neuro-jena.github.io/cat12-html/cat_starting.html)

Older version can be obtained [here](http://141.35.69.218/cat12/).

If you intend to install CAT12 from Github:
- Large files are stored with Git Large Storage, which you need to install first:
https://docs.github.com/de/repositories/working-with-files/managing-large-files/installing-git-large-file-storage
- Download ZIP CAT12:
https://github.com/ChristianGaser/cat12/archive/refs/heads/main.zip
- Remove any old CAT12 installation from spm12/toolbox
- Unzip cat12-main.zip to spm12/toolbox/ and rename it to spm12/toolbox/cat12

## Requirements
CAT12 is designed to work with SPM12 and Matlab versions 7.4 (R2007a) or newer. No additional toolboxes are required.

## Installation
- Remove the old cat12 folder in spm12/toolbox if existing
- Unpack the zip-file
- Copy the cat12 folder to the spm12/toolbox directory
- If once installed use the update function in CAT12 in order to check for new versions
- After restarting SPM12 either call CAT12 via the toolbox button or (as short-cut) type *cat12* on the Matlab command line.

## Download Standalone Version (no need for Matlab licene)
The advantage of the standalone version is that no Matlab license is needed. Only the (free) Matlab Runtime R2017b (v93) has to be downloaded. However, there are some limitations (e.g. no parallelization and no interactive help in the GUI version) and the standalone version is mainly intended to run without GUI on Unix systems. Please check the [ENIGMA CAT12](https://neuro-jena.github.io/enigma-cat12/#standalone) site fore more information and examples to call CAT12 from shell scripts.

The MATLAB Compiler Runtime (MCR) enables you to run applications compiled within MATLAB using MATLAB Compiler. MCR does not require a MATLAB license and can be used to run the MATLAB compiled program on computers which do not have MATLAB installed.

|CAT12 Standalone|MCR|
|---|---:|
<!--|[Mac](http://141.35.69.218/cat12/cat12_latest_R2017b_MCR_Mac.zip) |[Mac](https://ssd.mathworks.com/supportfiles/downloads/R2017b/deployment_files/R2017b/installers/maci64/MCR_R2017b_maci64_installer.dmg.zip)|
[Windows](http://141.35.69.218/cat12/cat12_latest_R2017b_MCR_Win.zip) |[Windows](https://ssd.mathworks.com/supportfiles/downloads/R2017b/deployment_files/R2017b/installers/win64/MCR_R2017b_win64_installer.exe)|
-->
[Linux](http://141.35.69.218/cat12/cat12_latest_R2017b_MCR_Linux.zip) |[Linux](https://ssd.mathworks.com/supportfiles/downloads/R2017b/deployment_files/R2017b/installers/glnxa64/MCR_R2017b_glnxa64_installer.zip)|

Please contact [me](mailto:christian.gaser@uni-jena.de) if you need versions for MacOS or Windows.