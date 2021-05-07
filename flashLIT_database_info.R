# flashLIT Database structure

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# flashLIT information

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# Causal inference informed, semi-automated systematic mapping for compex topics									

# Background ---									
# Systematic review challenged for complex topics (large numbers of interacting components, many papers, inconsistent terminology)									
# Automated grouping repeatable, but often not replicable (especially for adding new papers), and still requires interpretation of group themes and interactions.									
# Classic systematic review weighs sensitivity and specificity, narrowness of theme with volume of work. 									
# Term-list allocation problematic when incomplete list of terms known, lack of term sensitivity ot specificity for themes, lack of resources to manually resolve these.									
# More manual work = more chance for human error and bias.									
# Largely automated methods offer the benefits of improved systematic handling of literature for complex topics. These include revtools and scimeetr, which use different methods to cluster and summarize documents from a review. However, these provide different methods for 'natural' clustering, leaving the user to interpret these clusters, i.e. a descriptive clustering, rather than prescriptive. Further, these are repeatable, but not necessarily replicable with revised data.
# We propose a method that leverages the replicability of term-list allocation with automated grouping methods.								
# This can provide a living database that can be further queried and refined.									

# Aim	---								
# Develop a semi-automated process for systematic mapping of complex topics using a modern systematic review toolbox.									

# Process	---								
# 0a	Define the 'draft' causal inference map for the topic, at least defining main interacting themes as nodes.								
# 1a	Define 'library' using broad, catch-all search terms, compile using revtools								

# 2a	Automated grouping using revtools/scimeetr (inc references, authors). First scimeetr (auto group k) then use group k within revtools. Iteratively reduce to deliver groups with <10-20 members								
# 2b	If using more than one group method, compare group membership at each level between algorithms, identify 'robust' and 'uncertain' groups	
# 2c	(optional: Extract group-level keywords, characteristics, reading lists)								

# 3a	Create a 'node-term library' from nodes & known synonyms, expert knowledge, etc. Needs to be based on the causal paths/themes.								

# 4a	Create a term-document matrix for each of the title, keyword, abstract. These should contain terms! citations could be used in early iterations, but should not be needed as existing terms should lead to groups where other potential terms can be found. 								

# 5a	Map papers with enough support onto causal map nodes based on the node-term list (should start simple, but can develop into a specificity * occurence metric, for example).								
# 5b	For each node with new papers, investigate group-level keywords (2c) to supplement (3a), then repeat process. In more advanced repetitions, use automated methods to extract keyword-combinations								

# 6a	If the group has most members in the same node/theme, consider to classify othersin the group into that node/theme								
# 6b	At last resort, identify likely themes for each group manually via keywords (revtools, scimeetr), reading lists (scimeetr)								

# 7a	Classify likely node specificity and sensitiviy (likelihood of giving good information for that node) based on classification metadata
# 7b	Allocate likelihood for unallocated groups or papers.								

# 8a	Different versions can loop through parts of the process and be documented in metadata.	

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