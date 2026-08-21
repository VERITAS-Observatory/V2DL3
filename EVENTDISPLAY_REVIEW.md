# Review of the Eventdisplay path in V2DL3

Review date: 2026-08-04  
Reviewed revision: `448b2a3` (`eventdisplay-maintenance`, identical to `main` at review time)  
Scope: Eventdisplay conversion, shared FITS writers, index generation, comparison/query utilities, documentation, packaging, tests, and Eventdisplay CI. VEGAS-specific code and documentation were excluded.

Priority reassessment: 2026-08-21. The labels below reflect the production invocation that converts one explicitly paired anasum/effective-area file per run, selects exactly one response mode through `${m}`, always supplies a DB-FITS file and instrument epoch, explicitly selects the interpolator, writes to a run- and mode-specific output path, and may supply filtering/extrapolation options through `${V2DL3SELECT[@]}` and `${V2DL3OPT[@]}`.

The command shown does not invoke index generation, does not explicitly pass `--filename_to_obsid`, and uses `${m}` as one mode flag. Unless those options are hidden inside `${V2DL3OPT[@]}`, the related findings are not primary risks for this converter path. Conversely, defects in event selection, response interpolation, DB metadata, or FITS headers affect every generated run file and are prioritized accordingly.

## Executive summary

The normal Eventdisplay conversion path is compact and understandable: an anasum ROOT file and an effective-area ROOT file are converted into one DL3 FITS file, after which a separate command creates observation and HDU index tables. The implementation has useful separation between Eventdisplay ROOT extraction and shared FITS serialization, supports point-like and full-enclosure responses, and has a reproducible integration fixture.

The baseline is not yet ready for unattended scientific production without targeted fixes. The highest-priority risks for this command are unverified DB metadata applied to every run, event filters changing the run parameters used for IRF interpolation when `${V2DL3SELECT[@]}` contains `--evt_filter`, and exact-zero azimuth failure when the selected interpolator is the regular-grid implementation. Every output also carries an EDISP `CREF7` typo and a non-integer `MJDREFI` value that should be corrected or validated. The equal offset values (`THETA_LO == THETA_HI`) remain a low-priority compatibility question, not an established VERITAS response error. Filename-derived observation IDs, conflicting mode flags, and index-path handling are conditional or outside the shown converter invocation.

The repository's 57 unit tests pass, and newly generated baseline point-like and full-enclosure products match the stored references. Those facts establish reproducibility, not format or scientific correctness: the tests do not exercise most command-line behavior or semantic GADF validation.

## Typical documented workflow

The README describes the following Eventdisplay workflow:

1. Install with `pip install v2dl3`, or create `v2dl3Eventdisplay` from `environment-eventdisplay.yml` for development.
2. Convert one run with `v2dl3-eventdisplay -f RUN.anasum.root EFFECTIVE_AREA.root OUTPUT.fits`.
3. Select `--full-enclosure` when needed; point-like output is the default.
4. Optionally select an interpolator, allow extrapolation, clamp near an IRF boundary, add telescope multiplicity or database metadata, and filter events with YAML/JSON.
5. Run `v2dl3-generate-index-file` over the resulting per-observation files to create `obs-index.fits.gz` and `hdu-index.fits.gz` for Gammapy.

The implementation follows this broad flow:

- `pyV2DL3/script/v2dl3_for_Eventdisplay.py` defines the CLI.
- `EventDisplayDataSource` orchestrates event/GTI extraction and response generation.
- `eventdisplay/fillEVENTS.py` reads Eventdisplay trees, creates event metadata, and decodes the time mask.
- `IrfExtractor`, `IrfInterpolator`, and `eventdisplay/fillRESPONSE.py` select and interpolate IRFs.
- Shared `fillEVENTS.py`, `fillGTI.py`, and `fillRESPONSE.py` serialize FITS extensions.
- `generateObsHduIndex.py` and `script/generate_index_file.py` build the data-store indices.

## Findings

### ED-01 — Low: response offsets are serialized as sampled nodes

- [x] reassessed for the VERITAS/Eventdisplay production path. No established scientific error; retain as a compatibility check.

`find_camera_offsets()` extracts the unique `Woff` values from the effective-area file. For multiple offsets it returns the same array for both values passed to `THETA_LO` and `THETA_HI` (`eventdisplay/fillRESPONSE.py:143-164`). The response builders then evaluate the interpolator at each offset coordinate and store the resulting samples (`eventdisplay/fillRESPONSE.py:189-206`, `231-250`, and `332-342`). The local comment explicitly notes that these values may represent simulated points rather than interval boundaries.

The supplied VTS products do not show a V2DL3-specific defect:

- The generated VTS files use equal `THETA_LO` and `THETA_HI` values in their corresponding response extensions. These values match the `Woff` coordinates at which the Eventdisplay interpolator is evaluated.

Therefore, equal lower and upper values are being used here as a node representation, not necessarily as empty integration bins. The [GADF response specifications](https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/full_enclosure/aeff/) define these columns as the field-of-view offset axis but do not, by themselves, establish that every axis must have positive bin widths. Gammapy's current IRF serialization also writes identical LO/HI values for node-type axes ([source](https://github.com/gammapy/gammapy/blob/main/gammapy/irf/core.py)).

The original critical finding consequently conflated sampled nodes with interval edges. The comparison provides no evidence that the generated VTS offsets are wrong. A residual compatibility risk remains if a downstream consumer interprets these columns strictly as interval edges and attempts width- or solid-angle-based integration.

Recommended follow-up:

- Keep this as a low-priority format/consumer-compatibility check.
- Validate the generated files with the supported Gammapy reader and confirm that the offset coordinates and response values are loaded at the intended nodes.
- Do not convert the values to edges without first deciding the boundary policy and confirming that the target consumer requires interval semantics; doing so would change the meaning of the interpolated response grid.

### ED-18 — Medium and always-on: generated FITS metadata is inconsistent with its columns

The converter always writes an EDISP response HDU. Its `CREF7` keyword currently names the true-energy axis `ETRUE_LO:ETRUE_HI`, while the actual columns are `ENERG_LO` and `ENERG_HI` (`fillRESPONSE.py:45-49`). This is metadata-only, but it can confuse validators or readers that use `CREF7` to resolve the axis names.

The converter also always writes `MJDREFI` to EVENTS and GTI using `VTS_REFERENCE_MJD = 53402.0` (`constant.py:4`, `fillEVENTS.py:111-115`, and `fillGTI.py:26`). The integer part should be represented as an integer; the current header value is serialized as a floating-point value in generated products.

Suggested fix: change `CREF7` to `(ENERG_LO:ENERG_HI,MIGRA_LO:MIGRA_HI,THETA_LO:THETA_HI)` and define the reference MJD integer as `53402` (or cast explicitly at every header writer). Add a semantic FITS-header test to the normal point-like and full-enclosure CLI fixtures.

### ED-02 — Resolved: azimuth selection was not circular at north

- [x] fixed in #260 

Historical status: fixed in the current tree. `find_closest_az()` now normalizes the input azimuth and selects by circular angular distance; regression tests cover 359/0/360-degree and negative inputs. Retain the finding as closed unless production data show a regression.

### ED-03 — High when `--evt_filter` is used: event filtering can change the IRF applied to the observation

The event mask is applied before calculating mean pedestal variance, mean event altitude/azimuth, and the telescope mask (`eventdisplay/fillEVENTS.py:100-126`). Those masked values then define `NSBLEVEL`, `ALT_PNT`, `AZ_PNT`, and the zenith/pedestal/azimuth coordinates passed to IRF interpolation (`eventdisplay/fillEVENTS.py:53-68` and `73-85`).

This couples a user-level event selection to the instrument response. An energy, gammaness, or spatial filter can therefore change the IRF attached to the same observation. In the bundled sample, selecting 1–10 TeV changed the pedestal variance from about 7.2970 to 7.2952; the difference is small for that run but proves the coupling.

It also labels averages of reconstructed event directions as `ALT_PNT`/`AZ_PNT`, while RA/Dec pointing comes from `pointingDataReduced`. This is conceptually inconsistent and makes the header dependent on the event population.

Suggested fix:

- Derive IRF coordinates and pointing metadata from unfiltered run/pointing information.
- Apply `evt_filter` only to columns written to the EVENTS table.
- Calculate run-level pedestal variance from an unfiltered source, or document and validate the intended source explicitly.
- Test that adding a non-time event filter leaves all IRF arrays and observation-level headers unchanged.

### ED-04 — Deferred from the shown path: `--filename_to_obsid` creates invalid and inconsistent observation IDs

The option derives the ID with one `os.path.splitext()` call and writes it only to the EVENTS header (`v2dl3_for_Eventdisplay.py:160-166`). For `/tmp/custom-name.fits.gz`, this produces the string `custom-name.fits`, not an integer. The logging call formats that string with `%d`, producing a runtime logging traceback. The response extensions retain the original integer run number.

Reproduced output:

```text
EVENTS OBS_ID = 'custom-name.fits'
EFFECTIVE AREA OBS_ID = 64080
ENERGY DISPERSION OBS_ID = 64080
```

GADF defines the EVENTS `OBS_ID` as an integer unique observation identifier ([EVENTS specification](https://gamma-astro-data-formats.readthedocs.io/en/latest/events/events.html)). The index writer also declares `OBS_ID` as `>i8`, so arbitrary strings are incompatible with the next documented step.

Suggested fix: replace this option with an explicit integer `--obs-id`, validate it with Click, update the data source before generating any HDUs, and ensure every applicable extension and both indices use the same value. If filename inference is retained, strip all recognized FITS suffixes and require the remaining stem to parse as an integer.

### ED-05 — Deferred to index generation: HDU index paths are silently truncated

`gen_hdu_index()` fixes `FILE_DIR` to `S40` and `FILE_NAME` to `S54` (`generateObsHduIndex.py:63-78`). Astropy silently truncated a 60-character directory to 40 characters in a targeted reproduction. The resulting index can point to a nonexistent file while remaining a syntactically valid FITS table.

Suggested fix: determine widths from the maximum encoded values in the completed table, impose and validate any intentional format limit, and round-trip every indexed path in an integration test.

### ED-06 — Deferred from the shown path: both IRF mode flags are accepted and point-like silently wins

`--point-like` and `--full-enclosure` are independent flags (`v2dl3_for_Eventdisplay.py:31-40`). When both are true, response generation and serialization use `if point-like ... elif full-enclosure`, so the result contains only point-like IRFs (`eventdisplay/fillRESPONSE.py:375-447`; shared `fillRESPONSE.py:95-173`).

This was reproduced: the log reports both modes enabled, while the FITS file contains only point-like effective area and energy dispersion.

Suggested fix: make the modes a single `click.Choice` option or reject the conflicting combination with `click.UsageError`.

### ED-07 — Medium when `--evt_filter` is used: the documented energy-filter example fails

The README uses:

```yaml
ENERGY: [1, 10]
```

The filter is evaluated against raw Eventdisplay branch names before columns are renamed (`eventdisplay/fillEVENTS.py:259-282`). The actual branch is `Energy`, so the documented example raises `KeyError: 'ENERGY'`; `Energy: [1, 10]` works.

Suggested fix: define a stable public filter vocabulary matching output column names, translate it to ROOT branches internally, and validate unknown keys with a message listing supported names. At minimum, correct the README.

### ED-08 — Resolved in the current tree: valid filters that select zero events crash with an unrelated NumPy error

The current tree checks the selection result and raises `ZeroLengthEventList` before reductions such as `np.max` or `np.mean`. Keep the regression test because this is a direct failure mode when `${V2DL3SELECT[@]}` contains an event filter.

The historical failure was an unhelpful NumPy reduction error after a filter selected no events. The remaining filter-schema weakness is that malformed range lists and unknown keys can still fail with an unhelpful exception (`eventdisplay/fillEVENTS.py:269-274`).

Suggested follow-up: validate filter schema and keys before evaluation, and raise a domain-specific Click error identifying the filter and input file.

### ED-09 — Low/conditional: all-masked GTIs are handled, but missing time-mask fallback is fragile

The current `getGTI()` handles an all-zero mask and returns empty GTIs. The remaining issue is the missing-mask fallback: it catches `KeyError`, then immediately iterates through the same `timeMask` path (`eventdisplay/fillEVENTS.py:226-236`). If the whole `timeMask` node is absent rather than only `maskBits`, the fallback raises a second `KeyError` instead of using the run start/stop interval.

Suggested fix: make the fallback independent of the missing object and add a fixture with no `timeMask` node. Decide whether an all-bad run should produce an empty GTI or fail conversion with a clear quality message.

### ED-10 — Medium when `--force_extrapolation` is used: it does not provide the documented linear extrapolation with the default interpolator

The README states that this option linearly extrapolates the IRF. `check_parameter_range()` merely allows the out-of-range coordinate (`eventdisplay/fillRESPONSE.py:53-100`). With `RegularGridInterpolator`, `fill_value=None` does request linear extrapolation (`IrfInterpolator.py:135-145`). With the default `KNeighborsRegressor`, prediction remains a distance-weighted nearest-neighbor estimate (`IrfInterpolator.py:75-92`) and is not linear extrapolation.

This is a scientific semantics mismatch: the same CLI flag means different mathematics depending on another option.

Suggested fix: disallow forced extrapolation for KNN, implement and test an explicit extrapolator, or accurately document the method-specific behavior. Log the chosen algorithm and out-of-bounds distance in every extrapolated output.

### ED-11 — Resolved: zero-degree azimuth is valid in both interpolator paths

- [x] fixed in #261

`extract_irf()` now checks explicitly for `azimuth is None`, so `0.0` is accepted (`IrfExtractor.py:201-212`). `IrfInterpolator.interpolate()` now consistently requires the three coordinates used by all callers—`[pedvar, zenith, offset]`—without a special zero-azimuth branch (`IrfInterpolator.py:148-164`). Tests cover zero-degree extraction and regular-grid interpolation.

No further action is required for this finding.

### ED-12 — Deferred to index generation: index CLI existence checks ignore custom output names

When `--recreate` is absent, `generate_index_file.py:111-119` checks only hard-coded `obs-index.fits.gz` and `hdu-index.fits.gz`, not the values supplied through `--obs_index_file` and `--hdu_index_file`.

Consequences:

- Existing default-named files prevent creation even when different custom names were requested.
- Existing custom-named files are not detected and are overwritten because the lower-level writer always uses `overwrite=True`.

The converter itself also always overwrites its output (`v2dl3_for_Eventdisplay.py:167`) without an explicit `--overwrite` option.

Suggested fix: construct exact output `Path` objects once, check those paths, and require explicit overwrite authorization.

### ED-13 — Resolved: database metadata is now associated with and isolated from the run

- [x] fixed #263

The current implementation passes the anasum run number into `read_db_fits_file()` (`eventdisplay/fillEVENTS.py:71`). The DB reader now requires a supported run identifier column, selects exactly one matching DQM row, rejects missing/wrong/duplicate matches, exposes only supported auxiliary header fields, and rejects columns that would overwrite core Eventdisplay metadata (`eventdisplay/DBFitsFile.py`). `NSBLEVEL` and `QUALITY` remain derived from the anasum input, so DB metadata cannot alter the IRF query or run-quality value.

Regression tests cover matching-row selection, masked values, unsupported-column filtering, missing/ambiguous runs, core-metadata collisions, and propagation of the run number from the Eventdisplay event builder.

### ED-14 — Low for the CLI path: reusable Eventdisplay code depends implicitly on an active Click context

Range checks, interpolation setup, and extrapolation read CLI parameters through `click.get_current_context()` by default (`eventdisplay/fillRESPONSE.py:69-75` and `355-359`; `IrfInterpolator.py:135-139`). Calling the otherwise public-looking sequence `loadROOTFiles(..., "Eventdisplay"); datasource.fill_data()` outside the CLI reproduces `RuntimeError: There is no active click context`.

There is a hidden `use_click=False` path, but it is not a clean or documented API and option shapes differ (for example fuzzy-boundary handling).

Suggested fix: create a typed conversion configuration object in the CLI and pass it through the data-source/IRF layers. Domain code should not import Click or inspect global command context.

### ED-15 — Medium: production conversion relies on private Uproot internals and unchecked log text

Both `EventDisplayDataSource._fill_data_source_version()` and `eventdisplay_query_runparameters.py` access `members['fLines']._data`, where `_data` is private Uproot state. The version path checks that the log exists but not that the expected version line exists, so `[...][0]` can raise `IndexError` (`EventDisplayDataSource.py:87-95`). The query helper similarly assumes exact log phrases and split delimiters (`eventdisplay_query_runparameters.py:14-29`) and does not use a context manager.

Suggested fix: isolate ROOT log decoding behind one tested adapter, use public Uproot APIs where possible, handle absent/multiple matches explicitly, and include file/run context in errors.

### ED-16 — Resolved: comparison CLI now reports differences with a failing exit status

- [x] fixed in #264

`v2dl3-compareFitsFiles` still writes the FITSDiff report, but now exits with status 1 when `fd.identical` is false (`script/compareFitsFiles.py`). Regression tests cover identical and different FITS files. Tolerances and ignored keywords remain unchanged.

### ED-17 — Low: documentation and CLI text contain smaller inconsistencies

- [x] fixed in #265

- README says `--save_multiplicity` writes an `EVENTTYPE` keyword; code writes an `EVENT_TYPE` table column (`README.md:146`; shared `fillEVENTS.py:66-70`).
- The CLI help says “filter events form json or yaml” and the fuzzy-boundary help concatenates words because escaped source lines omit spaces (`v2dl3_for_Eventdisplay.py:64-87`).
- The README documents a tolerance argument as if singular, although the CLI expects repeated axis/value pairs.
- The index section links GADF v0.2 while files declare `HDUVERS=0.3` (`README.md:175-177`; `constant.py:12-15`).
- `environment-eventdisplay.yml` allows any Python `>=3.11`, while CI tests only 3.13 and package metadata claims Python 3.8 onward. The review environment resolved Python 3.14, which is not represented by the classifiers or CI.

Suggested fix: generate CLI examples from tested invocations, align supported Python constraints/classifiers/CI, and add a short Eventdisplay-specific reference page describing input tree expectations, response-mode semantics, filter names, output naming, and failure behavior.

## Cross-cutting design and quality recommendations

### 1. Add semantic DL3 validation, not only reference comparison

The current integration tests prove deterministic reproduction of stored output, including stored defects. Add a validator stage that independently checks:

- required extensions, columns, units, HDU class keywords, and integer observation IDs;
- consistent `OBS_ID`, time reference, telescope, and instrument metadata across extensions;
- finite axes with compatible shapes, and explicit semantics for interval axes versus sampled-node axes such as the Eventdisplay offset grid;
- non-negative effective area, normalized energy dispersion, and approximately normalized PSF;
- event times within observation/GTI bounds and monotonic ordering where required;
- successful loading and evaluation through the minimum and latest supported Gammapy versions;
- successful index round-trip from each `FILE_DIR`/`FILE_NAME` entry.

Keep reference comparisons as regression tests, but regenerate references only after semantic checks pass and require review of a summarized scientific diff.

### 2. Separate run metadata, event selection, and response construction

Create immutable models such as `ObservationMetadata`, `EventTable`, `GTI`, `IrfQuery`, and `ConversionOptions`. Extract observation metadata once from unfiltered run-level inputs, derive IRFs from that metadata, and apply event filters only to `EventTable`. This removes the current hidden data flow through mutable dictionaries and double-underscore attributes.

### 3. Validate at system boundaries

Use Click types/callbacks and explicit schema validation to reject:

- conflicting response modes;
- invalid/negative fuzzy tolerances;
- unsupported filter keys, malformed ranges, and empty selections;
- non-integer observation IDs;
- wrong-run or multi-row DQM inputs;
- IRF files lacking required branches, axes, or sufficient parameter coverage;
- output collisions unless overwrite is explicitly requested.

Errors should identify the file, run, branch/axis, supplied value, allowed range, and suggested correction.

### 4. Make interpolation behavior explicit and testable

Remove Click access from interpolation classes. Give each interpolation/extrapolation method a documented mathematical contract. Test exact grid nodes, interior points, boundaries, circular azimuth selection, single-dimension grids, empty bins, 0/360 degrees, fuzzy clamping, and out-of-range behavior. Record the method and query coordinates in reproducible provenance metadata or logs.

### 5. Strengthen the test matrix

Add unit and CLI tests for every finding above. Run them across the actually supported Python range. The Eventdisplay workflow currently ignores `README.md` changes in its pull-request path filter (`.github/workflows/v2dl3-eventdisplay.yml:5-15`), even though examples are part of the user interface; documentation changes should run at least CLI/example validation.

### 6. Make index generation transactional

Build and validate both tables in memory, write temporary files in the target directory, then atomically replace both final outputs only after both succeed. Validate duplicate observation IDs, metadata consistency, and path round-trips. Derive string widths from data rather than hard-coding them.

## Verification performed

- Read the Eventdisplay README workflow, actual `--help`, package/environment metadata, Eventdisplay modules, shared FITS writers, general utilities, tests, and Eventdisplay GitHub Actions workflow.
- Excluded all VEGAS implementation and VEGAS-only utilities from findings.
- Ran `pytest -q` in the `v2dl3Eventdisplay` environment: **59 passed**.
- Ran focused Flake8 fatal-error checks (`E9`, `F63`, `F7`, `F82`): no findings.
- Generated point-like and full-enclosure files from the bundled `64080` fixtures.
- Compared both outputs to stored reference FITS files: `FITSDiff` reported no differences.
- Ran Astropy FITS structural verification: no structural FITS errors.
- Independently checked response axis bounds and exposed the zero-width offset axes.
- Checked full-enclosure PSF values: finite, with numerical solid-angle integrals ranging roughly from 0.92 to 1.02 for the sample (useful as a baseline, but not a replacement for downstream-tool validation).
- Reproduced the azimuth wrap, filename observation-ID, README filter, no-Click-context, conflicting-mode, and fixed-width path issues described above; verified that the current tree handles zero-event selections and all-zero GTIs.

## Recommended implementation order

1. Correct EDISP `CREF7` and integer `MJDREFI`, then add semantic validation of the normal point-like and full-enclosure outputs.
2. If event filters are used in `${V2DL3SELECT[@]}`, decouple filtering from observation metadata and IRF coordinates; keep the zero-event regression test.
3. Fix exact-zero azimuth behavior for both interpolators and validate the explicitly selected `${INTERPOLATOR}` in batch runs.
4. If `${V2DL3OPT[@]}` enables extrapolation, correct KNN extrapolation semantics or disable that combination.
5. Harden the missing-time-mask fallback and document the current all-zero GTI policy.
6. Validate the sampled offset convention with the supported downstream reader; do not convert offsets to edges without a confirmed consumer requirement.
7. Defer filename-derived IDs, index path widths, custom-output checks, and transactional index writes unless those options/workflows are part of production.
8. Expand CLI/integration tests and correct the documentation.
