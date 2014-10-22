# HITChip atlas preprocessing

Custom preprocessing pipeline has been constructed for the HITChip atlas; to fetch and preprocess HITChip atlas on one go, use the following (check MySQL database (db) parameters; expect 30-60 minutes processing time). Only for internal use. Usage example:

<pre><code>
library(microbiome)
atlas <- FetchHITChipAtlas(allowed.projects = my.projects, 
		  dbuser = 'username', dbpwd = 'passu', dbname = 'Phyloarray', 
		  result.path = "results")

</pre></code>