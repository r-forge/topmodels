## Functions
### Forecasting functions (workhorse) 

Function name | S3 classes supported | S3 classes planned | TODOs
--- | --- | --- | ---
`procast()` | `lm`, `crch`, `disttree` | `gam`, `countreg`, `betareg`, `gamlss`, `bamlss` | many
`procast_setup()` | none | none | none
`newresponse()` | `default` | no | few
`qresiduals()` | `default` | no | many 

### Functions for graphical model assessment 


Class | Functions | `c()`, `rbind()` | `plot()` | `lines()` | `autoplot()` | TODOs
--- | --- | --- | --- | --- | --- | ---
`pithist` | `pithist.default()` | yes | yes | yes | yes | few
no | `qqrplot()` | no | no | no | no | many 
`reliagram` | `default`, `crch` | no | no | no | no | many
`rootogram` | `default`| yes | yes | yes | no | yes | few
`wormplot` | no | no | no | no | no | all


