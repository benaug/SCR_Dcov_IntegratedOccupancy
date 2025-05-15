# SCR_Dcov_IntegratedOccupancy

Spatial capture recapture model with spatial density covariates and a second observation model that is occupancy data (unmarked SCR data, too).
This is a closed population version of the model considered by Sun et al. (2019). There is a single and multisession version.

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2777

In this model, the occupancy data at each trap is related to the SCR process model following Ramsey et al. (2015).

https://wildlife.onlinelibrary.wiley.com/doi/full/10.1002/jwmg.851

Instead of probabilistically updating the latent capture histories for the occupancy data, I marginalize over individuals.
To speed up computation, I use the approach of Herliansyah et al. (2024) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

This model assumes occupancy data comes from same closed population as SCR data and detectors are not baited. 
For baited cameras, trap-level occupancy data can be used, but not trap by occasion occupancy data. 
I.e., did a trap detect anything or not over all occasions.
Also, it is assumed that SCR and occupancy detectors are not co-located where the same individual may be detected by both methods on the same site visit.

These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation