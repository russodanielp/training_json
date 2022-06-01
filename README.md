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
 
 2) Some responses are negative (e.g., -% inhibition, AID 411) which causes the hill equation to give all negatives 
 a. Could be as easy as just seeing of most of the concentrations are < 0 then multiplying by 1? 
 
 The Notebook `column_exploration.ipynb` displays what we're deadling with.  
 
 It appears the Box account has a file size limit of 15G, meaning the SQLite database
  would not fit.  
 
 ### Building database
 
 Okay, so, looks like the easiest path was just to parse all the JSON files for dose response information directly.  Because
  finding out which columns were "concentration" and which were "activity" through parsing of user submitted names, while 
  certainly doable, looks prone to information loss.
  
  `json_parsing.py`
  
  Instead, the `parse_json` json file attempts to parse the json file download from PuChemFTP into a format as returned
   by PugRest.  In the json files, some activities have a field `tc` which seems to identify the result as having a 
   concentration associated with it.  The `convert_json` function converts all the json files into the same (hopefully) 
   format as returned by PugRest (we could write some tests to confirm this).  The folder `json_from_csv` contains these
   csv files. 
   
   `curve_fitting.py`
   
   The function `build_sql_db()` loads the created csv files into a SQLite database.
   
   The function `add_concise_sql()` adds the corresponding concise csv files for these dose response assays.  
   This is important information to have.  They include active/inactive calls, ac50 values, etc.  
   
   `pug_rest.py`
   
   Two functions: `get_targets()` will get the target information for these assays and `load_targets_sql()` 
   loads them into a table in the SQLite database. 
   