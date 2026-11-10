# Transport scorecard

A separate scorecard for stratospheric tracer transport — age of air, residual circulation,
polar vortex and their variability — kept beside the two CESM tiers of
[EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md). It does not change any threshold or rule of
that protocol. It is the `jcm-strat` relative scorecard (phase 29), ported into this
repository.

---

## 1. Field contract

What a model has to supply. This extends the five fields of protocol §4.2 with the age tracer;
`Wzm` is optional here because no metric of this scorecard reads it.

| Name | Units | Definition | Metrics that read it |
|---|---|---|---|
| `Uzm` | m s⁻¹ | zonal-mean zonal wind | C1, D1–D4 |
| `Vzm` | m s⁻¹ | zonal-mean meridional wind | B1, B2, D4-w, D7 |
| `THzm` | K | zonal-mean potential temperature, reference 1000 hPa | B1, B2, D4-w, D7; C2 and D6 through T |
| `VTHzm` | K m s⁻¹ | zonal-mean **eddy** heat flux v'θ' | B1, B2, D4-w, D5, D6, D7 |
| `AOA150` | years | zonal-mean entry age: an age clock reset to zero wherever p > 150 hPa | A1–A4, D7 |
| `Wzm` | m s⁻¹ | log-pressure w̄ (protocol §4.2) — optional | none |
| `T` | K | zonal-mean temperature — optional; if absent, `T = THzm (p/1000)^κ` | C2, D6 |

- **Shape** `(time, plev, lat)`, daily means, every field on one common set of levels and
  latitudes.
- **`plev` in hPa**, any order; the reader sorts it ascending. Pressure must be known per level,
  independent of surface pressure, as in protocol §4.4.
- **`lat` in degrees**, global (−90° to +90°), any order; the reader sorts it ascending.
- **`time`** daily, each value labelled by the calendar day of its own mean. CF time with any
  calendar; a `noleap` calendar is read through `cftime`.
- **`VTHzm` is the eddy flux**, not the product of the zonal means — protocol §4.2 and
  convention 1 of [CLAUDE.md](../../CLAUDE.md). It must be formed from instantaneous fields
  before any time averaging; the jcm-strat runs form it within each 6-hourly snapshot and
  average four snapshots to a day.
- **`AOA150` is in years** (365.25 days). The jcm-strat files store the clock in days; their
  reader converts.
- **Levels**: the residual streamfunction Ψ* is integrated from the top supplied level down, so
  the levels have to reach well above 10 hPa (the jcm-strat runs reach 0.01 hPa, CESM h6
  4.5 × 10⁻⁶ hPa). The metrics read 5–150 hPa.
- **Record**: the window scored here is 2006–2010 (1826 days), because the references are
  matched to it (§3).
- Fill values must be masked before any reduction (protocol §5, rule 4).
