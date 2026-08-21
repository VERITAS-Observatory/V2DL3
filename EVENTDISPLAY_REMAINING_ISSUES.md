# Eventdisplay: remaining issues

Reviewed 2026-08-21 against `main...eventdisplay-maintenance` (merge base `79358e8`). This is a consolidated follow-up to `EVENTDISPLAY_REVIEW.md` and `CopilotEventdisplayReview.md`; issues already fixed on this branch are omitted.

Status after this implementation: all actionable issues below are addressed;
the `Reference-MJD type is wrong` item is intentionally excluded because it is
being fixed elsewhere.

## P1 — address before using the affected option in production

### Addressed — event selection changes the response for the same observation

`--evt_filter` is applied before the converter derives mean pedestal variance and mean reconstructed ALT/AZ. Those filtered values become `NSBLEVEL`, `ALT_PNT`/`AZ_PNT`, and the IRF query coordinates. Thus an energy, direction, or classifier cut can change the response attached to an otherwise identical run. The reconstructed-direction average is also labelled as pointing, while RA/Dec pointing is read from pointing data.

Implemented: derive run-level event metadata before applying the selection and use the selection only for `EVENTS`. Regression coverage verifies that the IRF query remains based on the full run.

References: `pyV2DL3/eventdisplay/fillEVENTS.py:40-87`, `pyV2DL3/eventdisplay/EventDisplayDataSource.py:31-40`.

### Addressed — `--filename_to_obsid` makes a non-integer, inconsistent observation ID

The option writes a filename-derived string only to the `EVENTS` header after all other HDUs have been built. Response HDUs keep the run number, and the index writer requires an integer `OBS_ID`. A `.fits.gz` filename also retains the `.fits` suffix after one `splitext` call.

Implemented without changing the option name: strip recognized FITS suffixes, require an integer stem, and propagate the resulting ID to every generated HDU.

Reference: `pyV2DL3/script/v2dl3_for_Eventdisplay.py:58-62,166-173`.

### Addressed — conflicting response-mode flags are silently accepted

Both `--point-like` and `--full-enclosure` can be supplied. Both flags are recorded as true, but response generation takes the point-like `if` branch, so the full-enclosure request is silently lost.

Implemented: retain both existing flags but raise a Click usage error when both are supplied.

References: `pyV2DL3/script/v2dl3_for_Eventdisplay.py:31-38,134-140`, `pyV2DL3/eventdisplay/fillRESPONSE.py:362-390`.

### Addressed — `--force_extrapolation` has two different meanings

The README promises linear extrapolation. That is true for `RegularGridInterpolator` (`fill_value=None`), but the default `KNeighborsRegressor` continues to produce a distance-weighted nearest-neighbour estimate. The same command-line option therefore has different scientific semantics depending on another option.

Implemented: reject forced extrapolation for KNN and document that the option requires `RegularGridInterpolator`.

References: `pyV2DL3/eventdisplay/IrfInterpolator.py:75-87,135-145`, `README.md:145-151`.

## P2 — correct for interoperable, reproducible products

### Reference-MJD type is wrong

`MJDREFI` is written from the floating constant `53402.0` to both `EVENTS` and `GTI`, although it declares the integer part of the reference MJD.

Fix: write `MJDREFI` as integer `53402` consistently. Add semantic header checks to point-like and full-enclosure integration tests.

References: `pyV2DL3/fillEVENTS.py:108-115`, `pyV2DL3/fillGTI.py:23-27`.

### Addressed — index generation can create bad or unexpectedly overwritten indices

The HDU index fixes `FILE_DIR` and `FILE_NAME` at 40 and 54 bytes, so longer paths are silently truncated. Separately, the CLI only checks the default index names before writing; requested custom names can be overwritten, while existing default names can block a custom request.

Implemented: size path columns from their encoded values, round-trip long paths in a test, and check the exact requested output paths before allowing overwrite.

References: `pyV2DL3/generateObsHduIndex.py:44-79`, `pyV2DL3/script/generate_index_file.py:111-119`.

### Addressed — Eventdisplay log parsing relies on private Uproot state and unchecked text

Version and run-parameter extraction access `fLines._data`, a private Uproot attribute, and assume an expected log line exists. A changed ROOT/Uproot representation or missing line can fail with an implementation exception (including an index error) rather than a contextual conversion error.

Implemented: centralize log decoding behind a tested adapter using public APIs where available, validate missing/ambiguous matches, and include input-file/run context in errors.

Reference: `pyV2DL3/eventdisplay/EventDisplayDataSource.py:78-95`.

## P3 — harden before broad reuse or unusual inputs

### Addressed — missing `timeMask` does not reach the intended fallback

The fallback catches a missing `maskBits`, then immediately dereferences the same `timeMask` node to list its keys. If `timeMask` itself is absent it raises another `KeyError` instead of using the run interval.

Implemented: make the fallback independent of `timeMask` and cover a missing-node fixture; a missing mask uses the full run interval, while an all-bad mask remains an empty-GTI result.

Reference: `pyV2DL3/eventdisplay/fillEVENTS.py:228-244`.

### Addressed — conversion internals require a live Click context

Response range checks and interpolator setup read the global Click context by default. The nominal data-source API consequently fails outside the CLI unless callers know the undocumented `use_click=False` escape hatch.

Implemented: non-CLI calls now default to explicit keyword options instead of reading a global Click context; the existing `use_click=True` compatibility path remains available.

References: `pyV2DL3/eventdisplay/fillRESPONSE.py:76-82,341-345`, `pyV2DL3/eventdisplay/IrfInterpolator.py:135-139`.

### Accepted — offset values are intentional sampled nodes

The output stores equal `THETA_LO` and `THETA_HI` values to represent sampled camera-offset nodes. This matches the current Eventdisplay interpolation and is covered by the response-array regression test. No conversion to interval edges is made because that would change the meaning of the sampled response grid.

The existing response-array test preserves these sampled coordinates for point-like and full-enclosure construction.

Reference: `pyV2DL3/eventdisplay/fillRESPONSE.py:350-360`.

## Verification

The full test suite passes: `71 passed`. `git diff --check` also passes.
