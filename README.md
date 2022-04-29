# PubChem AOP

Repository for code directed at AOP modeling using PubChem dose response assays

### Environment 

All code is written in Python and uses the `intro-chem` environment.  This can be 
created using the `environmental.yml` file. 


### Overview/Notes

1) PubChem assay contains many dose response assays, a list of which is available via PugRest 
(https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/type/doseresponse/aids/TXT).  The full output of this 
is stored in `data/dr_aids.txt`.  The dose response information from these can be retrieved via PugREST as well, 
but for very large assays (> 1,000 SIDs), you need to use a list key.  For that reason, there is a script `json_failed_aids.py`,
 that will parse the appropriate JSON file downloaded from PubChem's FTP and put it into csv format.  
 
 ### TODO
 
 1) Columns in the csv files are not standardized.  They can represent a diverse array of activity values (e.g., ac50 values, 
 hill slope values, or individual concentration responses).  None of these are standardized.  Need to:  
 a. Identify the columns that refer to a single concentration response.  
 b. Extract concentrations and units.    
 
 The Notebook `column_exploration.ipynb` displays what we're deadling with.  