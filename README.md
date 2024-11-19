# SCR_Dcov_IntegratedOccupancy

Spatial capture recapture model with spatial density covariates and a second observation model that is occupancy data.
This is a closed population version of the model considered by Sun et al. (2019). Assumes occupancy data comes from same
closed population as SCR data and detectors are not baited. For baited cameras, trap-level occupancy data can be used, but not trap by occasion occupancy data. I.e., did a trap detect anything or not over all occasions.
Also, it is assumed that SCR and occupancy detectors are not co-located where the same individual may be detected by both methods on the same site visit.
https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2777

Will post multisession versions when I get around to it.