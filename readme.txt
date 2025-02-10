"ppo.csv"
table consists of 216 observations of 14 variables where
observations are Myrmica scabrinodis workers sampled during the behavioural experiment
explanation of variables:
"id" - individual identification of every sampled worker, where the first part of the ID is the experimental colony code, the second is a worker number in this colony and the last part is the task performed by the individual - "F" forager, and "N" - intranidal worker
"pair" - number of a collected pair - each forager was collected with an intranidal worker 
"nest" - number of the experimental nest from which workers were sampled
"colony" - number of the mother colony from which workers came from
"task" - "N" - intranidal worker, "F" - forager
"rickia" - "1" - infected colony, "0" - uninfected colony
"lvl_inf" - average infection level of the colony
"age" - age of the worker when it was collected +- 2 days
"switch" - age when the forager was observed foraging for the first time
"plate" - number of spectrophotometric plate for phenoloxidase analysis
"po" - slope of the linear part of the spectrophotometric curve (Vmax) of the  samples without added chymotrypsin (active phenoloxidase)
"ppo" - slope of the linear part of the spectrophotometric curve (Vmax) of samples with added chymotrypsin (total phenoloxidase)
"head_len" - length of the worker's head
"head_wid" - width of the worker's head

"taskswitch.csv"
table consisting of 554 observations of 9 variables where
observations are all experimental workers
explanation of variables:
"mother" - number of the mother colony
"nest" - ID of the experimental nest
"nested" - number of the experimental nest from which workers were sampled
"task" - "nurse" - intranidal worker, 'forager'
"rickia" - "1" - infected colony, "0" - uninfected colony
"lvl_inf" - average infection level of the colony
"age" - age of the worker when it was collected +- 2 days
"switch" - age when the forager was observed foraging for the first time
"status" - 1 - switched task to forager, 0 - did not switch task to forager


"ppo_scaled.csv"
432 observations of 4 variables
observations are spectrophotometric measurements and variables are:
"age" - age of the worker when it was collected +- 2 days
"immune" - slope of the linear part of the spectrophotometric curve (Vmax)  
"type" - "po" - active phenoloxidase, 'ppo' - total phenoloxidase
"scaled" - for type = "ppo" divided by maxPO and multiplied by maxPPO, for type = 'po' left the same 
