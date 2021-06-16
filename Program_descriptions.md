# Description of the programs used in the MenA simulation project

## Utility files with functions for running the simulation and handling input, output, error-checking, and validation:

### ModelInputUtilities.R
  1. *GetMontaguDemogData:* given username, password, touchstone id, and destination folder, will download the 5 files used by this model from montagu.vaccineimpact.org using the API. (NOTE: authentication does not work from KPWA network computers. Tested from a personal laptop. To complete this step manually, download the following demographic files (names for touchstone &quot;201710gavi-5&quot;):
    1. 201710gavi-5\_dds-201710\_cbr\_both.csv
    2. 201710gavi-5\_dds-201710\_cdr\_both.csv
    3. 201710gavi-5\_dds-201710\_qq\_pop\_both.csv
    4. 201710gavi-5\_dds-201710\_tot\_pop\_both.csv
    5. 201710gavi-5\_dds-201710\_unwpp\_imr\_both.csv
  2. *GetDemographicParameters*: given the path to the downloaded files, a country and a date range, returns a dataframe of parameters for use in calling function. Note: option fillThreshold is intended to deal with the fact that some of the demographic files end at 2099, others at 2100, see function checkVIMCdates.
  3. *checkVIMCdates*: this function checks the date range specified against that in the downloaded files. If t exceeds the parameter date range, the function will fill the missing years with the closest existing year value, up to the number of years specified by fillThreshold
  4. *GetPopAgeDist*: given a country, start date, and directory of downloaded demographic files, this function calculated the fraction of the population in each off 7 age bands from the quinquennial year closest to the start date. It returns a vector with named elements.
  5. *GetVaccScenario*: given a country, scenario (routine or campaign) and the directory of downloaded files, returns a dataframe of yearly doses/coverage for use by simulation function. (NOTE: vaccine scenarios are not available via API at this time, so must be downloaded manually.)
  6. *GetDiseaseStateDist*: given region (hyperendemic or not) and directory containing file, returns a dataframe of disease state fractions by age group. File dist\_both.csv must be supplied in directory specified.
  7. *GetWAIFWmatrix*: given region (hyperendemic or not) and directory containing file, returns a dataframe of age-specific infection probabilities to be used by simulation function.
  8. *expandWAIFW* – input: waifw, expands to fit population age groups returns 4x361x2 matrix. The 3rd dimension is labelled rainy or dry so depending on the month, the simulation function will use the appropriate 3rd dimension to match up with the population matrix.
### MenA\_helper\_functions.R
  1. *GetInfectiveRatio* - input: population slice, output a vector of 4 ratios for 4 age groups?
  2. *Vaccinate* –input: given population slice, and vaccination dataframe (GetVaccScenario) , vaccination scenario, a date) function moves people into vaccinated category according to specifications if date is current.
### InitializePopulation.R
  1. One function: *Initialize population*. Initializes full age\*disease state\*time, populates the first time slice. Calls GetPopAgeDist and GetDiseaseStateDist to calculate starting values.
  2. INPUT: start &amp; end dates of simulation, popsize, country, region=&quot;not\_hyper&quot;
  3. OUTPUT: age\*disease state\*time matrix (361\*10\*duration in weeks) with first timepoint populated
### MenA\_OneSim.R
  1. One function: *MenASimulation*, Contains date loop of script described above, plus some setup
  2. INPUT: startdt, enddt, pop, fixedparams, countryparams, WAIFWmx, dxr
    1. Input population is the matrix, empty except for initial pop
    2. parameters are packaged into vectors and referred to by position
  3. RETURNS filled population x time matrix
  4. CALLS GetInfectiveRatio, Vaccinate
### MenA\_summarization\_functions.R
  1. *GetCohortSize*: sum total population by year of age (collapse 1-360, expand 361) INPUT: entire pop matrix (output of simulation), OUTPUT: dataframe [70 x years of simulation]
  2. *SummarizeOneSim*: sum cases by year of age, calculate deaths and dalys. INPUT: entire pop matrix (output of simulation), cfr, n (simulation counter, for labelling variable only.) OUTPUT: dataframe similar to above, with cases per age group then total deaths and DALYs
  3. *SummarizeForOutput* –INPUT : results list (outputs from SummarizeOneSm) , cohort (outputfrom getCohortSize), write (do you want to write a file? And if so, filename. Calculates mean cases, deaths, DALYs by year and age in years over all simulations,optionally writes csv file (directory is still hard coded.)
### MenA\_paramCheck.R
  1. *DemogNumVarExists*: Checks that the demographic inputs are correctly formatted.
  2. *IsCountryAndColAvailable*: Verifies that the requested country is in the data and the required columns are present.
  3. *GetFileName*: creates a name for an input file and verifies it is valid.
  4. *CheckDemogFileStructure*: Verifies demographic variables exist in specified file.
  5. *CheckDemogParameters*: Detailed checking of demographics
  6. *CheckSetParameters*: Make sure input parameters have valid values.

## Top-level program for running batches of simulations (MenA\_VaccSims.R)

1. Set up libraries, source function files
2. Set up variables to specify simulations to run (country, dates, vaccination scenario, etc.)
3. Load demographic parameters and vaccination data using functions in ModelInputUtilities.R
4. Package fixed parameters
5. Load WAIFW, using functions in ModelInputUtilities.R
6. Fill first slice of population matrix with starting values by calling function InitializePopulation.R
7. Make empty list for results
8. Execute simulation specified number of times
    1. Call MenASimulation
    2. Optional call to summarization functions GetCohortSize and SummarizeOneSim
1. After all simulations are complete
    1. Optional &quot;detail&quot; to write out results of first ten sims (for testing)
    2. Call SummarizeForOutput
