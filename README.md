## Changelog
### New features

#### Event Class Mode
- Using event class mode `-ec`, users may provide multiple effective areas per group. V2DL3 will infer event classes based on their effective area cut parameters and filter events to their respective event class based on those parameters.
- Outputs a FITS for each event class
- Easily extensible to define event classes based on more parameters, as well as multiple at once ([see dev notes](https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#Extending)).
- For now, this is only done based on the EA’s MSW ranges
- Users may use event class mode with a single EA as well to take advantage of the EA-based event filtering and cut validations provided by event class mode.
- Event class mode could be made default should these features be deemed appropriate for general use case.

#### PSF King Support
- V2DL3 now supports full-enclosure analysis by providing a king parameters file using `-k <king params file>` ([file format](https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#PSF_King_Parameters_File_Format))
- Outputs conform to psf_king format https://gamma-astro-data-formats.readthedocs.io/en/v0.1/irfs/psf/psf_king/index.html
- Note that psf king needs params to match your EA(s) MSW ranges

#### Spatial Exclusions
- The user may provide a file `-c <cutfile>` to define additional cuts for their stage6 run
- For now only supports spatial exclusion zones
- Easily extensible to more cuts using new `loadUserCuts()` in vegas/util.py ([see dev notes](https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#Extending_2))

#### FoVCuts
- Events can be automatically cut according to the EA's FoVCut parameters. 
- For now this only works in event class mode; will decoupled and added as a flag soon

#### Reconstruction Mode
- Users may choose ITM (.M3D) via `-r` flag (`-r 2` for M3D)
- Standard (.S) remains default if no choice is made

#### Cuts Validation
- Execution is stopped and the user is informed if problematic cut discrepancies are detected with the EA (e.g T^2 improper for point-like/FE) or between the EA and events (e.g events cut more deeply than an EA)
- Currently checks T^2 upper, MSW upper/lower, MaxHeight upper/lower, FoVCut upper/lower.
- May be overridden with `-o`
- Very easy to implement additional checks ([see dev notes](https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#Cuts_Validations))
- For now, only applies to event class mode to ensure that previous use cases are unaffected.

#### 4DMLM (`-4d`) Preset
- A convenient flag preset to take advantage of the new features for a 4D maximum likelihood method
- Presets `-l`, `-e`, `-k`, `-r 2`, `--full-enclosure` (runlist, event class mode, psf_king, ITM reco mode, and full enclosure)
- Takes `<runlist> <king params file>` as arguments

### Improvements

#### File pair mode bugfix
- File pair mode was broken for VEGAS due to event display kwarg `evt_filter`. VegasDataSource now supports receiving kwargs.

#### Improved exceptions and messages
- Solved `cpycppy` exception that would bury the user’s true exception under a mountain of cpycppy related messages.
- New exceptions to guard users and developers from using new features in undefined ways
- Added user and developer exceptions for more clear messages to the user in cases of program failure (e.g why AttributeError arises when using `–full-enclosure` without `psf_king`)

#### Improved Click support
- Recent versions of `click` cause program-halting exceptions that are difficult to solve. This patch supports both newer and older `click` releases.
