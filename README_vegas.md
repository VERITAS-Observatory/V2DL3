This README contains more detailed information about newer v2dl3-vegas options. See README.md for installation and basic usage.

## FoV Cutting

Events which fall outside the FoV defined in your EA file's cuts parameters will be automatically excluded from the resulting fits files. To disable FoV cutting, use `-nf`

## Reconstruction Mode

Use `-r 2` to choose ITM reconstruction instead of standard. If you're using a King PSF, make sure your parameters file is for the appropriate reconstruction type.

## Full Enclosure Analysis
`--full-enclosure`

v2dl3-vegas now supports full-enclosure via King PSF. Point-like remains default.

### King PSF
`-k <filepath>`

Some common parameter files are [available for download on the Wiki](https://veritas.sao.arizona.edu/wiki/V2dl3-vegas_dev_notes#King_PSF_Parameters_Files)

The axes values (azimuth, zenith, noise, etc) will be matched to the nearest PSF parameter line. So if you have runs near 40 zenith, for example, and your King PSF parameter values are mapped only for 0 and 20 zenith, you may want to get a new PSF parameters file made.

See the wiki link above for more information.

## Event Classes
`-ec`

Event Classes allow events of a run to be sorted to multiple fits files based on parameter(s).

When the user chooses `event_class_mode`, multiple effective area files may be provided for a group of runs, and these EA files will define multiple event classes based on their cuts parameters. [See Example](https://veritas.sao.arizona.edu/wiki/V2dl3-vegas_dev_notes#Event_Classes)

This is useful for splitting events from a run into multiple datasets for joint fitting. For example, disjoint signal and background classes based on MSW as developed by Alisha Chromey and Amanda Weinstein in "4DMLM" analysis [Example](https://github.com/VERITAS-Observatory/4DMLM-Analysis/blob/main/notebooks/crab-115-runs/JointCrab3DAnalysis2Events.ipynb)

Currently, event_class_mode only sorts to multiple event classes based on MSW intervals.
