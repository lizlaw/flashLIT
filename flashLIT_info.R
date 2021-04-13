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

# ENDSCRIPT