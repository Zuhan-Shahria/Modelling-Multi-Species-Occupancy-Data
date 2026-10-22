# Modelling Multi-Species Occupancy Data

**Effects of habitat type on the occupancy and detection of mammal communities in a logged Bornean forest**

A Bayesian multi-species occupancy analysis of camera trap data from the SAFE (Stability of Altered Forest Ecosystems) project in Sabah, Malaysian Borneo.

---

## Project Overview

Biodiversity is declining rapidly, and habitat loss from logging and agricultural expansion is a main driver. By the time a species is officially declared extinct, the chance to intervene has usually passed. Camera trap surveys combined with occupancy modelling offer a way to track which species use which habitats **before** populations reach the point of no return.

| | |
|---|---|
| **Data** | Camera traps at 260 sites, 2014-2017 |
| **Species** | 45 mammal species |
| **Habitats** | Trail, Clearing, Riparian, Forest, Disturbed Track |
| **Method** | Bayesian multi-species occupancy model (`msPGOcc`, `spOccupancy` in R) |
| **Inference** | Gibbs sampling (MCMC), 3 chains, 10,000 iterations each, 5,000 burn-in, thinning by 5 |

---

## The Questions

The central question:

> **How does habitat type affect the occupancy and detection of mammal species in a logged forest in Borneo?**

The challenge is that a camera failing to photograph a species does not mean the species is absent. So the real question is: *once we correct for imperfect detection, where do mammals actually live?*

This breaks down into:

1. **Occupancy:** Which habitats have the highest probability of being occupied by the mammal community as a whole?
2. **Species-level patterns:** Which individual species are associated with, or avoid, particular habitats?
3. **Detection:** How detectable is each species when it is present, and why do species differ?
4. **Method:** Does separating occupancy from detection change what we conclude compared with raw detection counts?

---

## Approach

The model separates two processes:

- **Ecological process (occupancy, ψ):** the probability that species *k* truly occupies site *i*, modelled with a logit link as a function of habitat type (Trail as the reference level).
- **Observation process (detection, p):** the probability a camera records species *k* in a given week *if it is present*, modelled as a species-specific constant.

Species are modelled jointly through hierarchical community-level priors, so rarely detected species borrow strength from the rest of the community. Priors are weakly informative (`N(0, 2.72)` on means, `Inverse-Gamma(0.1, 0.1)` on variances).

**Data preparation highlights**
- Kept only independent detection events (`New Contact`) to avoid burst-photo inflation.
- Matched detections to the correct deployment window, since some sites were deployed in multiple years.
- Split each deployment into weekly survey occasions.
- Consolidated 18 raw habitat labels into 5 ecological groups.
- Treated weeks with no active camera as `NA`, not as non-detections.

---

## What I Discovered

### 1. No single habitat dominates at the community level
Mean occupancy probability is roughly **0.25 in all five habitats**, with wide, heavily overlapping credible intervals. Opposing species-level responses cancel out when averaged, and Clearing and Disturbed Track had fewer deployments, adding uncertainty.

### 2. Some species stand out
- **Bearded Pig** is the most widespread species, with high occupancy in every habitat.
- **Malay Civet, Pig-tailed Macaque, Red Muntjac and Malay Porcupine** also show consistently high occupancy, which suggests they are habitat generalists.
- **Leopard Cat** shows a possible association with Forest, although its credible intervals overlap with other habitats, so this is suggestive rather than conclusive.
- **Banded Linsang, Sunda Clouded Leopard, Sunda Pangolin and several squirrels and treeshrews** have low estimated occupancy with wide intervals. The data cannot tell us whether they are truly rare or just hard to detect.

### 3. Detection is imperfect for every species
No species has a posterior mean weekly detection probability above **0.5**. Detection ranges from about 0.5 (Bearded Pig) to below 0.05 (Low's Squirrel, Sunda Clouded Leopard). Large, ground-dwelling species are easiest to detect; small arboreal species and cryptic carnivores are hardest. This is why raw detection counts would give a misleading picture and an occupancy model is needed.


---

## Decisions to Preserve Biodiversity

| # | Decision | Confidence | Why |
|---|---|---|---|
| 1 | **Protect a mosaic of habitats**, not one type | Moderate | No habitat is clearly best for the whole community; species use habitats differently |
| 2 | **Retain remaining forest fragments** in logged landscapes | Moderate (precautionary) | Leopard Cat and a few other species lean towards Forest; retention is low-cost insurance |
| 3 | **Fund targeted surveys** for data-poor species (Sunda Clouded Leopard, Marbled Cat, Banded Linsang, Sunda Pangolin) | High | Low detection means we cannot yet judge their status, and these are likely among the most threatened |
| 4 | **Deprioritise generalists** (Bearded Pig, Malay Civet, Pig-tailed Macaque) for habitat-specific protection | Moderate | They use all habitats widely; keep them as indicator species |
| 5 | **Follow up on Orangutan** before changing management | Moderate | Occupancy varies by habitat but the pattern is uncertain, and it is a flagship species |
| 6 | **Run repeated annual monitoring** analysed with occupancy models | High | The only way to catch declines early, before extinction becomes irreversible |
| 7 | **Do not change logging rules on this evidence alone** | Low | Logging roads and skid trails did not clearly reduce occupancy, but intervals are wide |

**One-sentence takeaway:** this landscape has no single habitat that protects every species, so keep habitat variety, retain remaining forest as a precaution, and invest in surveys for the rare species we currently cannot see.

---


## Tools

- R, `spOccupancy` (`msPGOcc`), `lubridate`, `ggplot2`
- Convergence checked with the Gelman-Rubin R-hat statistic


## References

- Doser, J.W., Finley, A.O., Kery, M. and Zipkin, E.F. (2022). *spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models.* Methods in Ecology and Evolution, 13(8), 1670-1678.
- MacKenzie, D.I. et al. (2002). *Estimating site occupancy rates when detection probabilities are less than one.* Ecology, 83(8), 2248-2255.
- Gelman, A. and Rubin, D.B. (1992). *Inference from iterative simulation using multiple sequences.* Statistical Science, 7(4), 457-472.
- WWF (2024). *Living Planet Report 2024.*
- IPBES (2019). *Global Assessment Report on Biodiversity and Ecosystem Services.*


## Author

**Zuhan Shahria**
