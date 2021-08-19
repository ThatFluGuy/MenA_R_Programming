# Simulating the impact of serogroup A meningococcal vaccines in Africa: COVID disruptions branch

_**Disclaimer**: This repository is under active development and may be incomplete._

## Contact
This repository was developed by Mike Jackson at the [Kaiser Permanente Washington Health Research Institute](https://www.kpwashingtonresearch.org/). For more infomation, please e-mail michael.l.jackson@kp.org.

## Background
Epidemics of bacterial meningitis caused by *Neisseria meningitidis* occur worldwide. The so-called "African Meningitis Belt" is a region of sub-Saharan Africa that covers parts of 26 nations, within which meningococcal epidemics occur with a frequency and intensity that is not seen elsewhere in the world. These epidemics are predominantly caused by serogroup A *N. meningitidis* (NmA), although epidemics can be caused by other serogroups as well. To combat these epidemics, the low-cost protein conjugate vaccine [MenAfriVac](https://www.path.org/articles/about-meningitis-vaccine-project/) (PsA-TT) has been developed for countries within the Belt. 

As the vaccine was first licensed, some key questions remained regarding vaccination strategy and target age groups. Our research team developed a mathematical model of NmA transmission and PsA-TT vaccination to address these questions.
```
Tartof S, Cohn A, Tarbangdo F, Djingarey MH, Messonnier N, Clark TA, et al.
"Identifying Optimal Vaccination Strategies for Serogroup A Neisseria meningitidis Conjugate Vaccine in the African Meningitis Belt"
PLoS One 2013; 8(5):e63605
```

This branch focuses on the potential impact of COVID-related disruptions on MenA vaccination programs. The purpose is to estimate the NmA cases/deaths that could result from decreased MenA vaccination coverage during COVID. This work has been published as:
```
Gaythorpe KAM, Abbas K, Huber J, Karachaliou A, Thakker N, Woodruff K, et al.
"Impact of COVID-19-related disruptions to measles, meningococcal A, and yellow fever vaccination in 10 countries"
eLife 2021;10:e67023
```


## Impact forecasting
Substantial funding for MenAfriVac, including mass campaigns and routine childhodd immunization, has come from [Gavi, the Vaccine Alliance](https://www.gavi.org/). As part of its mission, Gavi has worked to quantify the expected impact of the vaccination programs it funds. This effort has evolved into the [Vaccine Impact Modelling Consortium](https://www.vaccineimpact.org/). Our research team has participated in the VIMC since its inception, using our model to predict the expected number of NmA cases and deaths under various assumptions about vaccination coverage.
```
Li X, Mukandavire C, Cucunuba ZM, Londono SE, Abbas K, Clapham HE, et al.
"Estimating the health impact of vaccination against ten pathogens in 98 low-income and middle-income countries from 2000 to 2030: a modelling study"
The Lancet 2021; 387:398-408
```

As part of our work with VIMC we have refined our initial simulation model. The current version based on our most recent publication:
```
Jackson ML, Diallo AO, Medah I, Bicaba BW, Yameogo I, Koussoube D, et al.
"Initial validation of a simulation model for estimating the impact of serogroup A Neisseria meningitidis vaccination in the African meningitis belt"
PLoS One 2018; 13(10):e0206117
```

## Repository overview
This repository contains the R code and several ancillary files needed to run the epidemic simulation model. A description of the individual programs is provided in the file "Program_descriptions.md". Note that this repository does not contain all the data needed to run the simulations. Input files (including demographics and vaccine coverage) are stored in a separate directory that is not included in this repository. Those files include proprietary data on Gavi's demand forecasts, which cannot be publically posted here.
