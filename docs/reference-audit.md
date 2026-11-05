# Calculator bibliography review

Reviewed 2026-11-03 against catalog schemas from `a24694320f69b1cb0caabed5ffbbd7f82cace209` and source changes through `c1ec0fc13f979553a2fd399aacc4fef44557b4e8`.

## Result

- **123 retrieved BibTeX records:** 113 of the 125 distinct legacy citations, plus 10 background sources for DAG/Lean methods.
- **All 71 registered procedures** have at least one retrieved reference. Their schemas, calculation results, web reference lists and Word reports now use the imported records. The web result offers a BibTeX download for that procedure.
- **12 legacy citations remain unresolved.** They are excluded from the published bibliography and runtime reference lists, and retained in `references/unresolved.json` for follow-up. Absence of an export does not mean a work is fictitious.
- Google Scholar returned an unusual-traffic block. **None of these entries is represented as a Google Scholar export.** The sources are 119 Crossref content-negotiation exports, one Journal of Statistical Software export, one Springer export, and two existing institutional BibTeX records. The institutional exceptions are Lehmann (University of Geneva catalogue) and Hayes–Moulton (Johns Hopkins author record), retrieved through the public web index because direct access was denied.
- Citation identities were checked using titles, authors, editions/publication type and available journal details. This is a bibliography review, **not a validation of every formula, approximation, selectable method variant or numerical implementation**. Procedure-level references are not a claim that every variant has been separately verified.

## Source preservation

[`references.bib`](../references.bib) is assembled from [`references/exports/`](../references/exports/). No bibliographic fields were supplied from model knowledge or copied from the old citation strings. Citation keys are prefixed with stable inventory IDs to prevent collisions. Unneeded abstracts are omitted from the three exports containing them; all remaining fields are preserved. The two institutional records were extracted from existing BibTeX blocks, not reconstructed from prose. Their retrieval method is stated explicitly.

[`references/provenance.json`](../references/provenance.json) records each source URL, retrieval time, original key, retained-export SHA-256, and review notes. For downloaded exports with an omitted abstract it also records the original download hash. Display citations are formatted mechanically from those fields; HTML/TeX markup is removed only for display. Missing pages, source typos and electronic-publication dates are not silently replaced with the legacy text.

## Decisions that matter

| Record | Decision |
| --- | --- |
| Chow, Julious, Cohen, Ahn and Jennison–Turnbull books | Matched the intended editions; preserved source-reported publication dates even when they differ from a copyright year or the original print edition. |
| Diggle, Montgomery, Desu, Schlesselman and other books | Rejected chapter records, book reviews and wrong editions as substitutes for a whole-book citation. Some remain unresolved. |
| Blackwelder (1982) | Used the full article DOI `10.1016/0197-2456(82)90024-1`, not the similarly titled 1981 abstract. |
| Pocock et al., win ratio | Used European Heart Journal DOI `10.1093/eurheartj/ehr352`, not the conference abstract with a different author order. |
| Yu–Ganju, win-ratio sample size | Used original DOI `10.1002/sim.9297`, pages 950–963. A correction exists at `10.1002/sim.9388`; the implementation has not been checked against that correction in this bibliography task. |
| Hayes–Moulton, second edition | Rejected malformed Crossref BibTeX lacking title/authors; retrieved the existing BibTeX from the coauthor's institutional publication record. |
| McCullagh, ordinal regression | Preserved publisher pages 109–127; the old 109–142 range includes discussion. |
| Offen and Keene articles | Retrieved the full titles/subtitles and author lists rather than retaining abbreviated legacy text. |
| pmsampsize | Retrieved the official R-package DOI record. This does not also cite the separate Stata module. |
| Crossref journal dates and page gaps | Preserved first-online years and starting-page-only records. See per-record notes; these are limitations of the exports, not evidence that every old print-year/page-range citation was wrong. |
| Welch and Riley/lmerTest titles | Preserved source spelling/markup in BibTeX. Display formatting removes HTML, but does not alter the exported record. |

## Procedure coverage

Keys below resolve in `references.bib`. Unresolved IDs refer to the follow-up table. The mapping is maintained in `references/procedures.json`.

| Procedure | Retrieved keys | Unresolved legacy IDs |
| --- | --- | --- |
| `means.crossover` | `R001_Chow_2017`, `R002_Julious_2009`, `R003_Senn_2002` | — |
| `means.be.ratio` | `R004_Schuirmann_1987`, `R005_OWEN_1965`, `R001_Chow_2017` | R006 |
| `means.anova` | `R007_Cohen_2013`, `R009_Fleiss_1999`, `R010_Benton_2003` | R008 |
| `means.anova2` | `R007_Cohen_2013` | R011 |
| `means.one.t` | `R001_Chow_2017`, `R002_Julious_2009`, `R012_Machin_2018`, `R013_Lenth_1989`, `R010_Benton_2003` | — |
| `means.one.ni` | `R001_Chow_2017`, `R002_Julious_2009`, `R014_Blackwelder_1982` | — |
| `means.one.equiv` | `R001_Chow_2017`, `R002_Julious_2009`, `R004_Schuirmann_1987`, `R005_OWEN_1965` | — |
| `means.paired.t` | `R001_Chow_2017`, `R002_Julious_2009`, `R012_Machin_2018` | — |
| `means.two.repeated` | `R015_Diggle_2002`, `R016_Liu_1997` | — |
| `means.rm.anova` | `R007_Cohen_2013`, `R017_Muller_1989` | — |
| `means.two.slopes` | `R018_Fitzmaurice_2011`, `R015_Diggle_2002`, `R019_Ahn_2014` | — |
| `means.two.t` | `R001_Chow_2017`, `R002_Julious_2009`, `R012_Machin_2018`, `R013_Lenth_1989`, `R010_Benton_2003` | — |
| `means.two.welch` | `R001_Chow_2017`, `R002_Julious_2009`, `R012_Machin_2018`, `R013_Lenth_1989`, `R010_Benton_2003`, `R020_WELCH_1947` | — |
| `means.two.ni` | `R001_Chow_2017`, `R002_Julious_2009`, `R014_Blackwelder_1982` | — |
| `means.two.ancova` | `R021_Borm_2007`, `R001_Chow_2017` | — |
| `means.two.equiv` | `R001_Chow_2017`, `R002_Julious_2009`, `R004_Schuirmann_1987` | — |
| `means.coprimary` | `R022_Sozu_2011`, `R023_Sugimoto_2012`, `R024_Offen_2007` | — |
| `prop.cc` | `R026_Fleiss_2003` | R025, R027 |
| `prop.cc.matched` | `R028_Dupont_1988` | R025, R029 |
| `prop.chisq` | `R007_Cohen_2013`, `R010_Benton_2003`, `R030_Guenther_1977` | — |
| `prop.one` | `R001_Chow_2017`, `R026_Fleiss_2003`, `R012_Machin_2018`, `R031_Chernick_2002` | — |
| `prop.paired.mcnemar` | `R032_Connor_1987`, `R033_Lachin_1992`, `R026_Fleiss_2003` | — |
| `prop.trend` | `R034_Nam_1987`, `R035_Armitage_1955` | — |
| `prop.one.simon` | `R036_Simon_1989`, `R037_Jung_2004` | — |
| `prop.two` | `R026_Fleiss_2003`, `R001_Chow_2017`, `R012_Machin_2018`, `R038_Casagrande_1978`, `R039_Fisher_1935` | — |
| `prop.two.ni` | `R040_Farrington_1990`, `R041_Miettinen_1985`, `R001_Chow_2017`, `R002_Julious_2009` | — |
| `prop.two.equiv` | `R040_Farrington_1990`, `R041_Miettinen_1985`, `R001_Chow_2017` | — |
| `corr.one` | `R042_Guenther_1977`, `R007_Cohen_2013`, `R012_Machin_2018` | R043 |
| `corr.two` | `R007_Cohen_2013` | R043 |
| `reg.multiple` | `R007_Cohen_2013`, `R044_Gatsonis_1989`, `R010_Benton_2003` | — |
| `reg.slope` | `R045_Dupont_1998`, `R012_Machin_2018` | — |
| `reg.logistic` | `R046_Hsieh_1998`, `R047_Hsieh_1989` | — |
| `reg.prediction` | `R048_Riley_2020`, `R049_Riley_2018`, `R050_Riley_2018`, `R051_Riley_2020`, `R052_Ensor_2019`, `R053_Van_Houwelingen_1990` | — |
| `reg.cox` | `R054_Hsieh_2000`, `R055_Schoenfeld_1983` | — |
| `surv.logrank` | `R056_SCHOENFELD_1981`, `R057_Freedman_1982`, `R058_Lachin_1986`, `R012_Machin_2018` | — |
| `surv.logrank.ni` | `R056_SCHOENFELD_1981`, `R059_Jung_2005`, `R012_Machin_2018` | — |
| `rates.one.poisson` | `R060_GARWOOD_1936`, `R061_Gu_2008`, `R012_Machin_2018` | — |
| `rates.two.poisson` | `R061_Gu_2008`, `R012_Machin_2018` | — |
| `rates.two.negbin` | `R062_Zhu_2013`, `R063_Keene_2007` | — |
| `var.one` | `R012_Machin_2018` | R064 |
| `var.two` | `R012_Machin_2018` | R008 |
| `nonpar.signedrank` | `R065_Lehm:nonp:1975` | R066 |
| `nonpar.sign` | `R067_Noether_1987`, `R012_Machin_2018` | — |
| `nonpar.kruskal` | `R065_Lehm:nonp:1975`, `R068_Mahoney_1996` | — |
| `nonpar.mannwhitney` | `R067_Noether_1987`, `R065_Lehm:nonp:1975` | — |
| `gs.logrank` | `R069_Jennison_1999`, `R070_GORDON_LAN_1983`, `R071_O_Brien_1979`, `R072_POCOCK_1977`, `R073_KIM_1987`, `R074_Hwang_1990`, `R075_Wang_1987`, `R076_Pampallona_2001`, `R056_SCHOENFELD_1981` | — |
| `gs.means` | `R069_Jennison_1999`, `R070_GORDON_LAN_1983`, `R071_O_Brien_1979`, `R072_POCOCK_1977`, `R073_KIM_1987`, `R074_Hwang_1990`, `R075_Wang_1987`, `R076_Pampallona_2001` | — |
| `gs.props` | `R069_Jennison_1999`, `R070_GORDON_LAN_1983`, `R071_O_Brien_1979`, `R072_POCOCK_1977`, `R073_KIM_1987`, `R074_Hwang_1990`, `R075_Wang_1987`, `R076_Pampallona_2001` | — |
| `agree.icc` | `R077_Walter_1998`, `R078_Shrout_1979` | — |
| `agree.kappa` | `R079_Fleiss_1969`, `R080_Cantor_1996` | — |
| `diag.auc` | `R081_Hanley_1982`, `R082_Obuchowski_1998` | — |
| `diag.sens` | `R083_Flahault_2005`, `R084_Pepe_2003` | — |
| `diag.paired` | `R032_Connor_1987`, `R085_Alonzo_2002`, `R084_Pepe_2003` | — |
| `corr.ci` | `R012_Machin_2018` | — |
| `means.one.ci` | `R086_Kupper_1989`, `R012_Machin_2018` | — |
| `means.two.ci` | `R086_Kupper_1989`, `R012_Machin_2018` | — |
| `prop.one.ci` | `R087_Newcombe_1998`, `R088_CLOPPER_1934`, `R089_Agresti_1998`, `R026_Fleiss_2003` | — |
| `prop.two.ci` | `R090_Newcombe_1998`, `R091_Agresti_2000`, `R026_Fleiss_2003` | — |
| `means.two.cluster` | `R093_Campbell_2014`, `R012_Machin_2018` | R092 |
| `prop.two.cluster` | `R094_130b794513fd4698b61465bc34f695ac` | R092 |
| `cluster.sw` | `R095_Hussey_2007`, `R096_Hemming_2015` | — |
| `sim.lmm.cluster` | `R097_Eldridge_2006`, `R098_Kuznetsova_2017`, `R099_Landau_2012` | R092 |
| `sim.glmm.cluster` | `R100_Li_2015`, `R097_Eldridge_2006`, `R101_Pinheiro_1995`, `R099_Landau_2012` | — |
| `sim.sw` | `R095_Hussey_2007`, `R102_Hooper_2016`, `R103_Baio_2015`, `R099_Landau_2012` | — |
| `sim.winratio` | `R104_Pocock_2011`, `R105_Finkelstein_1999`, `R106_Bebu_2015`, `R107_Yu_2022`, `R099_Landau_2012` | — |
| `sim.lmm.longitudinal` | `R108_Laird_1982`, `R109_Satterthwaite_1946`, `R098_Kuznetsova_2017`, `R099_Landau_2012`, `R110_Arnold_2011` | — |
| `sim.mmrm` | `R111_Mallinckrodt_2001`, `R112_Lu_2008`, `R113_Anderson_1957`, `R099_Landau_2012` | — |
| `sim.means.robust` | `R067_Noether_1987`, `R114_Fagerland_2009`, `R020_WELCH_1947`, `R099_Landau_2012` | — |
| `sim.ordinal` | `R115_Whitehead_1993`, `R116_McCullagh_1980`, `R117_Roozenbeek_2011`, `R099_Landau_2012` | — |
| `sim.multiarm` | `R118_Dunnett_1955`, `R120_HOCHBERG_1988`, `R099_Landau_2012` | R119, R121 |
| `sim.surv.nph` | `R122_Fleming_2005`, `R123_Lin_2020`, `R124_Royston_2013`, `R056_SCHOENFELD_1981`, `R125_Genz_1992`, `R099_Landau_2012` | — |

## Unresolved records

These are the old strings under review, not retrieved BibTeX. They were not converted into new entries.

| ID | Existing citation | Why no entry was imported |
| --- | --- | --- |
| R006 | Diletti, E., Hauschke, D. & Steinijans, V. W. (1991). Sample size determination for bioequivalence assessment by means of confidence intervals. Int. J. Clin. Pharmacol. Ther. Toxicol. 29, 1–8. | No matching official BibTeX export found. The 1992 multiplicative-model paper is a different work. |
| R008 | Desu, M. M. & Raghavarao, D. (1990). Sample Size Methodology. Academic Press. | Google Books lists the book but its BibTeX download returned HTTP 403. Crossref candidates are book chapters. |
| R011 | Montgomery, D. C. (2017). Design and Analysis of Experiments, 9th ed. Wiley. | Google Books BibTeX download returned HTTP 403. Crossref matches were book reviews or different editions. |
| R025 | Schlesselman, J. J. (1982). Case-Control Studies: Design, Conduct, Analysis. Oxford University Press. | Google Books BibTeX download returned HTTP 403. Book-review DOIs were rejected. |
| R027 | Kelsey, J. L., Whittemore, A. S., Evans, A. S. & Thompson, W. D. (1996). Methods in Observational Epidemiology, 2nd ed. Oxford University Press. | Crossref book record uses the series title rather than the book title. No reliable complete export obtained. |
| R029 | Breslow, N. E. & Day, N. E. (1980). Statistical Methods in Cancer Research, Vol. I: The Analysis of Case-Control Studies. IARC. | Google Books BibTeX download returned HTTP 403. Book-review records were rejected. |
| R043 | Fisher, R. A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small sample. Metron 1, 3–32. | No matching official export found for the 1921 Metron paper. Related Pearson papers were rejected. |
| R064 | Ostle, B. & Malone, L. C. (1988). Statistics in Research, 4th ed. Iowa State University Press. | Crossref returned book reviews rather than the cited fourth edition. No matching export obtained. |
| R066 | Al-Sunduqchi, M. S. (1990). Determining the Appropriate Sample Size for Inferences Based on the Wilcoxon Statistics. PhD thesis, University of Wyoming. | Unpublished University of Wyoming dissertation. Citation appears in a publisher bibliography, but no official BibTeX record was retrieved. |
| R092 | Donner, A. & Klar, N. (2000). Design and Analysis of Cluster Randomization Trials in Health Research. Arnold. | Google Books BibTeX download returned HTTP 403; no exact Crossref book record retrieved. |
| R119 | Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6, 65–70. | No matching Crossref record retrieved. JSTOR citation access was blocked; related Holland–Copenhaver papers were rejected. |
| R121 | Horn, M. & Vollandt, R. (1995). Multiple Tests und Auswahlverfahren. Gustav Fischer (sample sizes for many-to-one comparisons). | Crossref matches are reviews of the book, not the book. No matching book export retrieved. |

## DAG and Lean scope

The ten additional records are catalogued in [`references/supplemental.json`](../references/supplemental.json) for g-formula, AIPW, TMLE, Super Learner, asymptotic inference and Kaplan–Meier background. They are included in the root bibliography, but **DAG and Lean result objects do not yet expose these references**. Their estimator/test/variance/nuisance-learner mappings still require separate implementation-level review.

## Validation

`python scripts/build_references.py --check` parses every retained export, verifies its hash and DOI (where available), rejects duplicate keys/DOIs, checks every mapping, and compares the generated BibTeX and C++ data against their committed versions. C++ tests check reference coverage for every registered procedure and matching schema/result exports. Neither check certifies the statistical appropriateness of a citation.
