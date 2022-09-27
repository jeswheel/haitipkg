# haitipkg

## TODO: 

- [ ] Always check/update function documentation. 
- [ ] Fix package dependencies:
  - [ ] Make sure that the package depends on `pomp`, but possibly also `panelPomp`, and `spatPomp`.
- [ ] Make sure error of installing package from `GitHub` using `devtools::install_github()` is fixed (in particular, this error was primarily when trying to view the function documentation in `R-studio`. 
- [ ] When installing the package, there is the following warning that should be addressed: 
``` 
Warning: replacing previous import ‘pomp::mcap’ by ‘spatPomp::mcap’ when loading ‘haitipkg’
```
I don't think anyone is using either of these functions, but it should be looked at more closely so that other users don't get this warning.
- [ ] Some parameters missing from documentation of the function `covars`

### Model 1:


### Model 2: 


### Model 3:

### Complete: 

- [x] Refine and test `project_from_filter2`, and delete `project_from_filter`.
- [x] Refine `fit_haiti2()` and `est_logLik2()` as needed.
- [x] Create function `fit_haiti1()` that will fit Model 1 from scratch. 
- [x] Document the datasets in the package.
- [x] Make project_from_filtering compatable with spatPomp/panelPomp models.
- [x] Create TODO list
- [x] Create function `fit_haiti2()` that will fit Model 3 from scratch. 
- [x] Create function `est_logLik2()` that will estimate the log-likelihood for model 2.
- [x] Create function `est_logLik3()` that will estimate the log-likelihood for model 3.
- [x] Create function `fit_haiti3()` that will fit Model 3 from scratch. 

