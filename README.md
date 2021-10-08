# R package MortCast

[![R Build Status](https://github.com/PPgp/MortCast/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/PPgp/MortCast/actions?workflow=R-CMD-check)

### Estimation and Projection of Age-Specific Mortality Rates

Age-specific mortality rates are estimated and projected using 
the Kannisto, Lee-Carter and related methods as described in [Sevcikova et al. (2016)](https://link.springer.com/chapter/10.1007%2F978-3-319-26603-9_15).

The main functions are:

* **cokannisto**: Extrapolates given mortality rates into higher
          ages using the Coherent Kannisto method. The original
          Kannisto method (with sex-independent extrapolation) is
          avalable in the function **kannisto**.
* **lileecarter.estimate**: Estimates the coherent Lee-Carter
          parameters for male and female mortality rates (Li and Lee
          2005), i.e. sex-independent parameters a<sub>x</sub> and k<sub>t</sub>, and the
          coherent parameter b<sub>x</sub>. In addition, it computes the ultimate
          b<sub>x</sub><sup>u</sup> for rotation (Li et al. 2013).  The underlying
          sex-independent estimation is implemented in the function
          **leecarter.estimate**.
* **mortcast**: Using estimated coherent Lee-Carter parameters
          and given future sex-specific life expectancies, it projects
          age-specific mortality rates, while (by default) rotating the
          b<sub>x</sub> parameter as described in Li et al. (2013).

Functions contained in the package can be used for both, 5-year and 1-year age groups.

Other methods for forecasting mortality rates are available:

* **pmd**: pattern of mortality decline

* **mlt**: model life tables, using [UN lookup tables]({https://www.un.org/development/desa/pd/data/extended-model-life-tables}). Note that the tables were updated in version 2.6-0 (October 2021). For previous version of the tables install 2.5-0, e.g. `devtools::install_github("PPgp/MortCast@v2.5-0")`

* **logquad**: log-quadratic mortality model

* **mortcast.blend**: combining two different methods

A life table can be constructed using the **life.table** function.

#### References 

Li, N. and Lee, R. D. (2005). Coherent mortality forecasts for a
     group of populations: An extension of the Lee-Carter method.
     Demography, 42, 575-594.

Li, N., Lee, R. D. and Gerland, P. (2013). Extending the
     Lee-Carter method to model the rotation of age patterns of
     mortality decline for long-term projections. Demography, 50,
     2037-2051.
