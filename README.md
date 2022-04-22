# haitipkg

## TODO: 

- [ ] Make [R/get_elimProbs.R](R/get_elimProbs.R) to work for model 2.
- [ ] Make project_from_filtering compatable with spatPomp/panelPomp models.
- [ ] Always check/update function documentation. 
- [ ] Fix package dependencies:
  - [ ] Make sure that the package depends on `pomp`, but possibly also `panelPomp`, and `spatPomp`.
- [ ] Make sure error of installing package from `GitHub` using `devtools::install_github()` is fixed (in particular, this error was primarily when trying to view the function documentation in `R-studio`. 
- [ ] When installing the package, there is the following warning that should be addressed: 
``` 
Warning: replacing previous import ‘pomp::mcap’ by ‘spatPomp::mcap’ when loading ‘haitipkg’
```
I don't think anyone is using either of these functions, but it should be looked at more closely so that other users don't get this warning.
- [ ] Document the datasets in the package.
- [ ] Some parameters missing from documentation of the function `covars`

### Model 1:

- [ ] Create function `fit_haiti1()` that will fit Model 1 from scratch. 
- [ ] Create function `est_logLik1()` that will estimate the log-likelihood for model 1.

### Model 2: 

- [ ] Refine `fit_haiti2()` and `est_logLik2()` as needed.

### Model 3:

- [ ] Refine and test `project_from_filter2`, and delete `project_from_filter`.
- [ ] Combine `fit_haiti{i}()` and `est_logLik{i}()` for each $i$ into single functions `fit_haiti()` and `est_logLik()`, respectively. 



### Complete: 

- [x] Create TODO list
- [x] Create function `fit_haiti2()` that will fit Model 3 from scratch. 
- [x] Create function `est_logLik2()` that will estimate the log-likelihood for model 2.
- [x] Create function `est_logLik3()` that will estimate the log-likelihood for model 3.
- [x] Create function `fit_haiti3()` that will fit Model 3 from scratch. 
