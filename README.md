# <em>InvaCost</em> Damage Cost:Management Expenditure Ratio
<a href="https://doi.org/10.5281/zenodo.10801171"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10801171.svg" alt="DOI"></a>
<img align="right" src="www/InvaCostLogoIdea10.jpg" alt="insect damage icon" width="150" style="margin-top: 20px">

National-level assessement of the drivers of invasive-species costs (ratio of damage costs:management expenditure & rate of ratio change) based on socio-economic traits of countries (cost data derived from the <a href="https://github.com/Farewe/invacost"><em>InvaCost</em></a> database)

Coding error discovered in 2026 and now updated (see <a href="https://github.com/cjabradshaw/InvaCostDamMgmRatio/tree/main/scripts">Scripts</a> folder)

<br>
Prof <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
April 2023; updated January 2024 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>
contributors: <a href="https://www.researchgate.net/profile/Philip-Hulme-2">Phil Hulme</a>, <a href="https://carleton.ca/biology/people/emma-hudgins/">Emma Hudgins</a>, <a href="https://www.mcgill.ca/qls/researchers/brian-leung">Brian Leung</a>, <a href="https://portal.findresearcher.sdu.dk/en/persons/mkour">Melina Kourantidou</a>, <a href="https://www.cee-m.fr/member/courtois-pierre/">Pierre Courtois</a>, <a href="https://scholar.google.com/citations?user=59VAYs4AAAAJ&hl=en">Anna Turbelin</a>, <a href="https://www.trinity.edu/directory/smcdermo">Shana McDermott</a>, <a href="https://www.uidaho.edu/cals/agricultural-economics-and-rural-sociology/our-people/katherine-lee">Katie Lee</a>, <a href="https://www.linkedin.com/in/danish-ali-ahmed-655934192/">Danish Ahmed</a>, <a href="https://www.research.ed.ac.uk/en/persons/guillaume-latombe">Guillaume Latombe</a>, <a href="https://azimpremjiuniversity.edu.in/people/alok-bang">Alok Bang</a>, <a href="https://www.abdn.ac.uk/people/thomas.bodey/">Thomas Bodey</a>, <a href="https://scholar.google.com/citations?user=fwHUGm0AAAAJ&hl=de">Phillip Haubrock</a>, <a href="https://www.flinders.edu.au/people/frederik.saltre">Frédérik Saltré</a>, <a href="https://www.ese.universite-paris-saclay.fr/en/team-members/franck-courchamp/">Franck Courchamp</a><br>
<br>
Accompanies paper:<br>
<br>
Bradshaw, CJA, PE Hulme, EJ Hudgins, B Leung, M Kourantidou, P Courtois, AJ Turbelin, SM McDermott, K Lee, DA Ahmed, G Latombe, A Bang, TW Bodey, PJ Haubrock, F Saltré, F Courchamp. 2024. <a href="http://doi.org/10.1016/j.ecolecon.2024.108166">Damage costs from invasive species exceed management expenditure in nations experiencing lower economic activity</a>. <em>Ecological Economics</em> doi:10.1016/j.ecolecon.2024.108166

## Abstract
While data on biological invasions and their economic toll are increasingly available, drivers of susceptibility to damage and cost-effectiveness of management in reducing long-term costs remain poorly understood. We used data describing the damage costs of, and management expenditure on, invasive species among 56 nations between 2000 and 2020 reported in the <em>InvaCost</em> database to test the overarching hypothesis that higher-income nations and those with higher trade volume have a higher efficiency to limit the damage incurred by invasive species by spending relatively more on management. We also tested whether nations with (<em>i</em>) more corruption have a reduced capacity to manage invasive species, leading to relatively higher damage costs, (<em>ii</em>) more educated citizens or greater technological and scientific output allow for improved incentives and ability to manage invasive species, thereby reducing relative damage costs, and (<em>iii</em>) economies based on higher primary resource dependencies (e.g., agriculture) are at greater risk of incurring high costs of invasive species, and so all other conditions being equal, have higher relative damage costs compared to management expenditure. By focusing on the ratio between damage costs and management expenditure, we analyse the willingness of countries to invest in management as a function of the extent of the damage suffered. We show that economic activity, measured by the volume of trade, is the main determinant of this ratio — the greater the volume, the smaller the ratio. We also found a higher rate of increase in the damage:management ratio as a country’'s proportion of total land area devoted to agriculture increased, suggesting that a higher economic dependency on agriculture predisposes a country to greater damage costs from invasive species over time. When considering the proportion of total costs identified as damage-related, results indicated that higher government investment in education produced higher proportional damage, and lower corruption and lower trade volume both reduced proportional damage. Our overall results suggest that wealthier nations with high per-capita imports of goods and services are more susceptible to damage, but also have a greater capacity to reduce it, and are therefore less threatened by biological invasions than countries with fewer resources and lower imports.
<br>
<br>
Based on (now out-of-date) preprint:<br>
<br>
Bradshaw, CJA, PE Hulme, EJ Hudgins, B Leung, M Kourantidou, P Courtois, AJ Turbelin, SM McDermott, K Lee, DA Ahmed, G Latombe, A Bang, TW Bodey, PJ Haubrock, F Saltré, F Courchamp. <a href="http://doi.org/10.2139/ssrn.4587717">Weaker economies experience higher relative damage costs arising from biological invasions</a>. <em></em>SSRN</em> doi:10.2139/ssrn.4587717
<br>

## <a href="https://github.com/cjabradshaw/InvaCostDamMgmRatio/tree/main/scripts">Scripts</a>
- <code>ntlScaleCostDiffInvaCostGithub2024.R</code>: main R code for analysis
- <code>new_lmer_AIC_tables3.R</code>: source code for information-theoretic algorithms
- <code>r.squared.R</code>: source code for calculating goodness-of-fit for linear models (including mixed-effects models)

## <a href="https://github.com/cjabradshaw/InvaCostDamMgmRatio/tree/main/data">Data</a>
- <em>GDPpc.csv</em>: per capita gross domestic product by country (source: <a href="https://data.worldbank.org/indicator/NY.GDP.PCAP.CD">World Bank</a>)
- <em>CPI.csv</em>: corruption perception index (source: <a href="https://www.transparency.org/en/cpi/2021">Transparency International</a>)
- <em>govexpedu.csv</em>: government expenditure on all education (% of GDP; source: <a href="https://data.worldbank.org/indicator/SE.XPD.TOTL.GD.ZS">World Bank</a>)
- <em>faoag.csv</em>: value added proportion of GDP from agriculture, fisheries, and forestry (source: <a href="https://www.fao.org/faostat/en/#data/MK">Food and Agriculture Organization</a> of the United Nations)
- <em>GHSI2022.csv</em>: global health security index (source: <a href="https://www.ghsindex.org/report-model/">Global Health Security Index</a>)
- <em>pcAgrLand.csv</em>: % land surface area devoted to agriculture (source: <a href="https://data.worldbank.org/indicator/AG.LND.AGRI.ZS">World Bank</a>)
- <em>importGS.csv</em>: imports of goods and services (source: <a href="https://data.worldbank.org/indicator/NE.IMP.GNFS.CD">World Bank</a>)
- <em>stjarticles.csv</em>: scientific & technical journal articles (source: <a href="https://data.worldbank.org/indicator/IP.JRN.ARTC.SC">World Bank</a>)
- <em>pop2021.csv</em>: 2021 national population size (source: <a href="https://data.worldbank.org/indicator/SP.POP.TOTL">World Bank</a>)
- <em>continent.countryINVACOST.csv</em> & <em>fao.cntry.code.csv</em> are code files for merging countries and regions across different datasets


## Required R packages
<code><a href="https://cran.r-project.org/web/packages/invacost/invacost.pdf">invacost</a></code>, <code>lme4</code>, <code>dismo</code>, <code>gbm</code>, <code>boot</code>, <code>VIM</code>, <code>mice</code>, <code>performance</code>, <code>sjPlot</code>, <code>rworldmap</code>, <code>rgeos</code>, <code>SpatialEpi</code>, <code>nlme</code>, <code>rcompanion</code>, <code>jsonlite</code>

## Corrected temporal reanalysis

The original temporal-resampling predicate in
`scripts/ntlScaleCostDiffInvaCostGithub2024.R` incorrectly used `>=` for
both interval bounds. It now uses the stated lower and upper bounds. The
resulting historical script retains the original interval endpoints, which
share boundary years.

`scripts/corrected_temporal_reanalysis.R` is the reproducible post-publication
reanalysis. It expands eligible observed, high-reliability costs to annual
records, uses a fixed seed, handles empty country-period strata explicitly,
propagates bootstrap uncertainty from paired damage and management resamples,
uses ML for information-criterion comparisons, and evaluates the pre-specified
combined model set across eight imputations rather than averaging imputed
predictors. It writes ignored outputs, including input and package metadata:

```sh
Rscript scripts/corrected_temporal_reanalysis.R outputs/corrected-temporal-reanalysis 1000
```

It evaluates inclusive bounded three-year windows (for a minimal comparison
with the historical endpoints), disjoint three-year windows
(2000--2002, ..., 2018--2020), and relaxed disjoint five- and seven-year
windows. Wider windows improve within-window coverage at the cost of temporal
resolution, so they are sensitivity analyses rather than replacements for the
stated three-year method. The generated `run_metadata.csv` records the exact
package and data version.

Using `invacost` 1.1.7 (whose bundled database is InvaCost v4.1), the
three-year definitions retained 48--50 countries for the cross-sectional
ratio, but only 23--26 for the temporal rate. Their lowest mean-AICc model was
intercept-only. Relaxing the windows changed the preferred ratio model:
five-year bins selected agricultural land, whereas seven-year bins selected
imports (pooled inclusion weight 0.74). This sensitivity to the temporal
definition means that the relaxed bins cannot provide robust support for the
published imports conclusion; nor do they remedy the sparse temporal-rate
sample.

## External validation data and analyses

`data/un_comtrade_pathway_imports_2016_2020.csv` is a committed UN Comtrade
snapshot of annual all-import, live-animal, fish/aquatic-product, live-plant
and wood-product imports for the analysis countries. It covers 46 countries
with complete or partial 2016--2020 data. `data/griis_country_species_richness.csv`
is a GBIF snapshot of GRIIS country-checklist richness, retaining the GBIF
dataset UUID for each of the 47 matched countries. Refresh either snapshot
with:

```sh
Rscript scripts/download_external_snapshots.R
```

`scripts/external_validation_analysis.R` tests the pathway and total-import
predictors, their interaction with the existing corruption-capacity covariate,
GRIIS richness, separate log damage and management outcomes, and a
source-reference reporting-effort sensitivity:

```sh
Rscript scripts/external_validation_analysis.R
```

For the seven-year sensitivity outcome, pathway models use 42 common complete
cases. Wood-product imports had the lowest AICc for the ratio and damage
models (weights 0.32 and 0.34), while all imports ranked second for the ratio
(weight 0.16). Management conditional on damage had no clearly dominant
model. Crucially, the ratio model including all imports and reporting effort
had weight 0.79, compared with 0.07 for imports alone. These exploratory
results are consistent with an ascertainment-sensitive association, not a
robust, stand-alone total-imports effect.

## Eurostat environmental-capacity subset

`data/eurostat_environmental_protection_expenditure_2016_2020.csv` is a
committed snapshot from Eurostat table `gov_10a_exp`: general-government
COFOG GF05 environmental-protection expenditure as a percentage of GDP. It is
explicitly a broad environmental-capacity proxy, **not** invasive-species
management expenditure. Refresh the snapshot and rerun the expanded
validation analysis with:

```sh
Rscript scripts/download_eurostat_environment_proxy.R
Rscript scripts/external_validation_analysis.R
```

Twelve countries report at least one 2016--2020 value, and ten overlap the
pathway-trade ratio subset. The intercept-only model had AICc weight 0.77;
the environmental-capacity and all-imports models had weights 0.12 and 0.10,
respectively. This small, geographically restricted proxy analysis therefore
does not corroborate the proposed imports--capacity mechanism and must not be
treated as validation of national invasive-species management expenditure.

## Adaptive boosted regression tree sensitivity

`scripts/adaptive_brt_external_validation.R` applies the adaptive BRT
calibration strategy used in the
[`wealthwellbeingageingpop`](https://github.com/cjabradshaw/wealthwellbeingageingpop)
analysis. For the external-validation ratio dataset it searches bag fractions,
learning rates, tree complexities, fold counts and admissible node sizes, then
records every calibration attempt and bootstrap failure rather than treating
non-convergence as a result. It assesses total imports, wood-product imports,
corruption-based capacity, GRIIS richness and reporting effort on the 42
complete cases:

```sh
Rscript scripts/adaptive_brt_external_validation.R 30
```

The 30-replicate run produced 29 successful, acceptable bootstrap fits. The
calibrated model had cross-validated correlation 0.41. Median relative
influences were 31.4% for GRIIS richness, 20.2% for wood-product imports,
19.7% for reporting effort, 14.5% for capacity and 13.3% for total imports;
all 95% bootstrap intervals were broad. This is exploratory, observational
evidence and not confirmation of a robust total-imports effect.

<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="200" style="margin-top: 20px"></a>
<a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="200" style="margin-top: 20px"></a>
