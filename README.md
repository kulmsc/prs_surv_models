# prs_surv_models

This repository contains all of the scripts that go along with a manuscript entited: 
"Benchmarking Polygenic Risk Score Model Assumptions: towards more accurate risk assessment"

The key idea behind the manuscript is that many of the models employed to assess polygenic risk scores contain flaws that can be directly mitigated with known modifications.
We split these flaws into five catagories, each of which corresponded to a group of model variations that we hypothesized would perform better (or at least differently from) a plain model.
Each of the model variations were used to predict eigteen different diseases, each of which had a corresponding polygenic risk score calculated from the PGS Catalog.

Our findings show that flawed models can lead to significantly improper estimates of a polygenic risk score's accuracy.  While the nature of these flaws makes it difficult to determine if the change in accuracy better or worse represents reality, we make a strong case that the model variations are improvements (and provide ample supporting references).

These scripts can be referenced to enact your own model variations.
