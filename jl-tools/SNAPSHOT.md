jl-tools snapshot from commit 9d776fb, 2026-06-08

## Local modifications (after snapshot)

The following edits were made to the vendored modules *after* copying from the source
repo, to support the deposited reproduction scripts. They diverge from the snapshot
commit above and should be upstreamed before the next re-snapshot.

- **abaqus.jl** — `get_connectivity_from_inp`: changed the file-path construction from
  `join([dir, fn])` to `joinpath(dir, fn)`, so callers can pass directories without a
  trailing separator. Behavior-preserving for existing callers.
- **makiePolyPlot.jl** — `plotNodalValues`: added optional `fig_descriptor` fields,
  each with a fallback so existing callers are unaffected:
  `clim` (explicit color limits; falls back to symmetric/auto), `z_aspect` (Axis3 z
  aspect ratio; default 0.5), `lift` (z-offset for `node_coords_plot` markers; default 0).
- **faultOperations.jl** — `tractions_from_sij_and_n`: node-id column (col 1) is now
  written as Integer rather than Float (cosmetic; read path is unaffected).


### Revision of the paper (November 2026)
- **faultOperations.jl** — `tractions_from_sij_and_n`: optional `dir_path` and `fname`
  arguments typed `Union{String,Nothing}` (previously `String` with a `nothing` default,
  which raised a `MethodError` when omitted); `reduce_matrix_with_repeated_ids` called
  with `out_matrix_type="Int-Float"` so that non-integer stresses are averaged correctly.
  Added `check_t0_tractions` (sanity check at t0: compressive normal tractions and
  |τ0| ≤ μ_max |σ′n0|, default μ_max = 0.61), called from the traction scripts.
- **makiePolyPlot.jl** — `plotNodalValues`: added further optional `fig_descriptor`
  fields, each with a fallback so existing callers are unaffected: `equal_xy` (equal
  horizontal scale; default false), `scalebar_km`/`scalebar_xy` (scale bar in data
  coordinates), `north_arrow`/`north_xy` (north arrow), `fontsize_annot` (annotation and
  colorbar font size; default 18), `coast_linewidth` (default 1), and `clip_lines`
  (clip `line_data` to the plotted extent; default true).