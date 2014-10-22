## Database computer

To access HITChip database from WUR database computer, do the following:

1. Login to the new database machine. It will ask two passwords during login. You can get these from the [admins](contact)

1. Create a new output folder for profiling files. To do this, open the folder icon on the left panel, navigate to the "Users" directory and Create New Folder with right-mouse click.

1. Open terminal window (the black screen icon on the left panel).

1. Start R by giving the standard command 'R' in the terminal window.


### HITChip

To extract HITChip data, see [these instructions](hitchip).


### MITChip

To extract MITChip data, give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "mit", dbpwd = "passu", dbname = "phyloarray_mit")
```

This prompts you through profiling options described in [full
instructions](protocol).


### PITChip (older chip)

To extract PITChip data (older chip), give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "pit", dbpwd = "passu", dbname = "phyloarray_pit")
```

This prompts you through profiling options described in [full
instructions](protocol).



### PITChip (new chip)

To extract PITChip data (new chip), give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "pit", dbpwd = "passu", dbname = "pitchipdb")
```

This prompts you through profiling options described in [full
instructions](protocol).


### ChickChip (older chip)

To extract ChickChip data (older chip), give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "mit", dbpwd = "passu", dbname = "chickchipdb")
```

This prompts you through profiling options described in [full
instructions](protocol).




