# nirstorm
Brainstorm plugin for fNIRS data analysis. 

Current features include classical within-subject analysis comprising motion-correction, MBLL and window-averaging.
## Authors

 * ''Thomas Vincent, PERFORM Centre and physics dpt., Concordia University, Montreal, Canada''
 * ''Zhengchen Cai, PERFORM Centre and physics dpt., Concordia University, Montreal, Canada''
 * ''Alexis Machado, Multimodal Functional Imaging Lab., Biomedical Engineering Dpt, McGill University, Montreal, Canada''
 * ''Louis Bherer, Centre de recherche EPIC, Institut de Cardiologie de Montréal, Montréal, Canada''
 * ''Jean-Marc Lina, Electrical Engineering Dpt, Ecole de Technologie Supérieure, Montréal, Canada''
 * ''Christophe Grova, PERFORM Centre and physics dpt., Concordia University, Montreal, Canada''

## Installation

[Brainstorm](http://neuroimage.usc.edu/brainstorm/) must be installed prior to installing nirstorm. It can also help to have brainstorm running while installing nirstorm.

The script `nst_install.m` takes care of copying or linking processes and functions into the brainstorm user folder.

Parts of the nirstorm plugin may already be shipped with the lastest brainstorm version and are available in the process selection menu in the "NIRS" submenu (modified bear lambert law and bad channel tagging).
The current installation will override them.

### Copy installation (windows, linux)

To copy all processes and functions into the brainstorm user folder, run under matlab:
```matlab
>> nst_install('copy');
```
When updates are downloaded, this installation command must be run again for changes to take effect.

### Linked installation (linux only)

To create symbolic links of all processes and functions into the brainstorm user folder, run under matlab:
```matlab
>> nst_install('link');
```
When updates are downloaded, this installation command has to be run again only if there are new files.

## Usage

The main documentation is in the form tutorials available on the [nirstorm github project wiki](https://github.com/Nirstorm/nirstorm/wiki#tutorials).
