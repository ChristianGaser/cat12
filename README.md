# CAT: Computational Anatomy Toolbox
This toolbox is an extension to [SPM](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Wellcome Department of Cognitive Neurology) to provide computational anatomy. This covers diverse morphometric methods such as voxel-based morphometry (VBM), surface-based morphometry (SBM), deformation-based morphometry (DBM), and region- or label-based morphometry (RBM).

It is developed by Christian Gaser and Robert Dahnke (Jena University Hospital, Departments of Psychiatry and Neurology) and free but copyright software, distributed under the terms of the [GNU General Public License](http://www.gnu.org/licenses/gpl-2.0.html) as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

## Download
[CAT toolbox](http://dbm.neuro.uni-jena.de/cat12/cat12_latest.zip)

[CAT manual](https://neuro-jena.github.io/cat12-html)
[Getting Started](https://neuro-jena.github.io/cat12-help/#get_started)
[Quick Start Guide](https://neuro-jena.github.io/cat12-html/cat_starting.html)

Older version can be obtained [here](http://dbm.neuro.uni-jena.de/cat12/).

If you intend to install CAT from Github:
- Download ZIP CAT:
https://github.com/ChristianGaser/cat12/archive/refs/heads/main.zip
- Remove any old CAT installation from spm/toolbox
- Unzip cat12-main.zip to spm/toolbox/ and rename it to spm12/toolbox/CAT

## Requirements
CAT is designed to work with SPM and Matlab versions 7.4 (R2007a) or newer. No additional toolboxes are required.

## Installation
- Remove the old cat12 or CAT folder in spm/toolbox if existing
- Unpack the zip-file
- Copy the CAT folder to the spm/toolbox directory
- If once installed use the update function in CAT in order to check for new versions
- After restarting SPM either call CAT via the toolbox button or (as short-cut) type *cat12* on the Matlab command line.

## Download Standalone Version (no need for Matlab licene)
The advantage of the standalone version is that no Matlab license is needed. Only the (free) Matlab Runtime R2023b (v232) has to be downloaded. However, there are some limitations (e.g. no parallelization and no interactive help in the GUI version) and the standalone version is mainly intended to run without GUI on Unix systems. Please check the [ENIGMA CAT](https://neuro-jena.github.io/enigma-cat12/#standalone) site fore more information and examples to call CAT from shell scripts.

The MATLAB Compiler Runtime (MCR) enables you to run applications compiled within MATLAB using MATLAB Compiler. MCR does not require a MATLAB license and can be used to run the MATLAB compiled program on computers which do not have MATLAB installed.

|CAT Standalone|MCR|
|---|---:|
[Linux](https://github.com/ChristianGaser/cat12/releases/latest/download/CAT_R2023b_MCR_Linux.zip) |[Linux](https://ssd.mathworks.com/supportfiles/downloads/R2023b/Release/10/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2023b_Update_10_glnxa64.zip)|
[Mac Intel](https://github.com/ChristianGaser/cat12/releases/latest/download/CAT_R2023b_MCR_Mac.zip) |[Mac](https://ssd.mathworks.com/supportfiles/downloads/R2023b/Release/10/deployment_files/installer/complete/maci64/MATLAB_Runtime_R2023b_Update_10_maci64.dmg.zip)|
[Mac ARM64](https://github.com/ChristianGaser/cat12/releases/latest/download/CAT_R2023b_MCR_Mac_arm64.zip) |[Mac ARM64](https://ssd.mathworks.com/supportfiles/downloads/R2023b/Release/10/deployment_files/installer/complete/maca64/MATLAB_Runtime_R2023b_Update_10_maca64.dmg)|
[Windows](https://github.com/ChristianGaser/cat12/releases/latest/download/CAT_R2023b_MCR_Win.zip) |[Windows](https://ssd.mathworks.com/supportfiles/downloads/R2023b/Release/10/deployment_files/installer/complete/win64/MATLAB_Runtime_R2023b_Update_10_win64.zip)|

<!--Please contact [me](mailto:christian.gaser@uni-jena.de) if you need other versions.-->