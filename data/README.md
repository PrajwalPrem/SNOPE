# Data

Put your S2 astrometry + radial-velocity tables here (or point
`data_path_pos` / `data_path_rv` at wherever they live).

Expected format -- both are **headerless** CSVs:

`tab_gillessen_pos.csv`
```
t_yr, alpha_mas, alpha_err_mas, delta_mas, delta_err_mas
```

`tab_gillessen_vr.csv`
```
t_yr, v_los_kms, v_los_err_kms
```

These are not redistributed here -- use your own copy of the
Gillessen et al. (or GRAVITY Collaboration) S2 astrometry/RV tables,
or any other compatible dataset.
