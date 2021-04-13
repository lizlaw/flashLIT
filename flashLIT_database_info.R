# flashLIT Database structure

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com


# Database structure---			

# Keys:									
# PID 	Paper ID								
# GID	Group ID								
# TID	Term ID								
# NID	Node ID		Nodes must be nested within themes, but nodes can == themes. in simplest form, 1 node is the whole library.						
# Theme	From theme list		Themes should = paths through the causal graph. 						
# RID	Researcher ID								
# AID 	Additional term ID. Additional terms are e.g. review, covariates, etc.							

# Paper details---									
# PID	Authors	Title	KW	Abstract	References	Year	Journal	Citation	DOI    (from revtools/import, need to make tools use same input)
# PID	doi exists,	date,	doi2text,	n sections,	methods extractable			(DOI check)
# PID	method	GID, level		(grouping)					
# PID	node	node certainty	node method	nodecheck method	RID	source	date		
# PID	theme	theme certainty	theme method	themecheck method	RID	source	date		

# Term document matrices---									
# PID - TID	n	(term document matrix)							
# PID - AID	n	(additional term document matrix to potentially identify e.g. lab, field, experiment, review, metaanalysis, etc.)							

# Group characteristics ---									
# GID	n	keywords	some sort of metric indicating overlap between grouping methods						
# GID	theme	RID	method	(groups identified to theme - if done)					
# GID	node	RID	method	(groups identified to node - if done)		
# GID readinglist type, PIDs allocated as this (readinglists e.g. of review, highly cited)

# Term-node identifiers	---
# theme node
# category list (for Aterms)
# TID	term	node	node specificity	RID	method	source	date	
# AID	Aterm	category	these can be identifying e.g. review papers, or species, locations/habitats etc. Can initiate with e.g. country lists, common habitat types, common bee groups and species. 						
# 
# Node and theme characteristics	--- (created from the above, used for queries)								
# NID	number of papers, reading lists (scimeetr), number of likely review papers.								
# Theme	number of papers, reading lists (scimeetr), number of likely review papers.								

# Version metadata ---									
# VID	Description	included PID, GID, TID, AID, Theme, NID, and resulting PID:NID lists. 							

# ! There are likely to be several versions of many of these tables. Therefore dates and ways to sort and split them are required, including a 'discarded' version of them.									
# I suggest there be a default (pre-compiled) version									

# Classification model for missing PID: nodes can be based on the term-document matrix and/or group membership.									

# Ways to utilize the database ---								
# Search for papers on a specific topic, using nodes and theme information (and combinations of nodes)									
# From a known paper, identify others that are likely similar (group, keywords, terms, nodes, themes, cite, cited by)									
# Identify the number of papers attributed to each node and various node-combinations (arcs)									
# Identify nodes/arcs(i.e. node pairs)/theme knowledge (number of papers) and number of 'review' papers									
# Identify key knowledge gaps  (note importance of these will depend on bayesian network, but perhaps potential to get estimate based on graph theory).									
# Identify potential for meta-analysis for well-informed nodes (consistent, specific use of terms, number of papers, covariates, existing meta-analyses and reviews)									
# Note current doi2text coverage (doi, html access, and check of import quality, for future potential use in extracting terms from methods.									
# Methodological details to experiment with ---
# Group methods: How consistent are the groups between methods, in different themes, and with different year sets? 
# Group similarity can be: for each paper i, how many papers in Gxi are also in Gyi (group x method and group y method).
# Term-node-classifier: this is the hinge in allocating the papers to nodes. I propose by starting only with very specific node terms, even if these are broader themes. Summary could be sum or max. Threshold could be a tuning parameter. Likely more than one mention required, at least in abstract, but mention in title likely suitable by itself. 
# Group methods - search locations: include authors/citations in defining groups or not? Authors likely specialize. Citations should have many more relevant refs, but also open can of worms. 
# Term-lookup - include citations or not? I think initially not. If including will need to determine a higher thresholds, likely these depend on topics, cultures
# ** many of these might require good documentation on where terms come from. 
# Summaries can be as a network (all paths between nodes) or as pathways (allowing for all pairwise arcs along the path).

# Develop as an R package. --- (would be ideal, but for now, series of functions)
# Tools to create
# Tools to query
# Database (including term lists)

# How to have user interaction? How to have users contribute papers, and suggested nodes, and changes? ---
# Suggest main version should be 'starting data' and 'updates' available for each of the terms
# Updates can then use version number lists + updates to replicate data.

# ENDSCRIPT