#import "common.typ": *

= Introduction <introduction>

Although most #ndbd interactions in the CUORE detector are contained within a single crystal, approximately 5% deposit energy across adjacent crystals. This analysis isolates these dual-site events to constrain the #ndbd decay signature. Using the CUORE Background Model (BM) analysis with the 1 ton-yr data release as the starting point, we aim to extract a constraint on the #ndbd half-life leveraging a Bayesian Markov Chain Monte Carlo optimization method. Restricting the analytical phase space to a limited two-dimensional energy range suppresses coincident backgrounds, yielding competitive relative sensitivity compared to the primary single-site constraint.

This analysis differs from the standard CUORE signal extraction in that it uses a two-dimensional fit over a thin diagonal band in $(E_1, E_2)$, where the available data are comparatively sparse. It also compares this restricted region to a background model constrained over a much wider energy range, so residual differences in spectral shape are more visible (Ex #Co60 sum peak), and because the fit is built from simulated spectra rather than a fully analytic model, nuisance effects cannot be incorporated in a purely closed-form way.

@data and @fit_model detail data selection and simulation inputs. @methods outlines analysis methods including the Bayesian Markov Chain Monte Carlo (MCMC) methodology used to derive the #ndbd half-life posterior distribution. @results and @systematics present the fit results and systematic uncertainties.
