# InvaCost Damage:Management Cost Ratio
<img align="right" src="www/InvaCostLogoIdea10.jpg" alt="insect damage icon" width="150" style="margin-top: 20px">

National-level assessement of the drivers of invasive-species costs (ratio of damage:management costs & rate of ratio change) based on socio-economic traits of countries (cost data derived from the <a href="https://github.com/Farewe/invacost"><em>InvaCost</em></a> database)

<br>
Prof <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
April 2023 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>
contributors: <a href="https://www.researchgate.net/profile/Philip-Hulme-2">Phil Hulme</a>, <a href="https://carleton.ca/biology/people/emma-hudgins/">Emma Hudgins</a>, <a href="https://www.mcgill.ca/qls/researchers/brian-leung">Brian Leung</a>, <a href="https://www.cee-m.fr/member/courtois-pierre/">Pierre Courtois</a>, <a href="https://scholar.google.com/citations?user=59VAYs4AAAAJ&hl=en">Anna Turbelin</a>, <a href="https://www.trinity.edu/directory/smcdermo">Shana McDermott</a>, <a href="https://www.uidaho.edu/cals/agricultural-economics-and-rural-sociology/our-people/katherine-lee">Katie Lee</a>, <a href="https://www.linkedin.com/in/danish-ali-ahmed-655934192/">Danish Ahmed</a>, <a href="https://azimpremjiuniversity.edu.in/people/alok-bang">Alok Bang</a>, <a href="[https://www.abdn.ac.uk/people/thomas.bodey/](https://www.research.ed.ac.uk/en/persons/guillaume-latombe)">Guillaume Latombe</a>, <a href="https://www.abdn.ac.uk/people/thomas.bodey/">Thomas Bodey</a>, <a href="https://scholar.google.com/citations?user=fwHUGm0AAAAJ&hl=de">Phillip Haubrock</a>, <a href="https://www.ese.universite-paris-saclay.fr/en/team-members/franck-courchamp/">Franck Courchamp</a><br>
<br>
Accompanies paper:<br>
<br>
Bradshaw, CJA, PE Hulme, EJ Hudgins, B Leung, M Kourantidou, P Courtois, AJ Turbelin, SM McDermott, K Lee, DA Ahmed, G Latombe, A Bang, TW Bodey, PJ Haubrock, F Saltré, F Courchamp. Weaker economies experience higher relative damage costs arising from biological invasions (in review)
<br>
<br>

## <a href="https://github.com/cjabradshaw/InvaCostDamMgmRatio/tree/main/scripts">Scripts</a>
- <code>ntlScaleCostDiffInvaCostgithub.R</code>: main R code for analysis
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
- <code>invacost</code>
- <code>lme4</code>
- <code>dismo</code>
- <code>gbm</code>
- <code>boot</code>
- <code>VIM</code>
- <code>mice</code>
- <code>performance</code>
- <code>sjPlot</code>
- <code>rworldmap</code>
- <code>rgeos</code>
- <code>SpatialEpi</code>
- <code>nlme</code>
- <code>rcompanion</code>

<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="200" style="margin-top: 20px"></a>
<a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="200" style="margin-top: 20px"></a>
