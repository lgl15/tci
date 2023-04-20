# TCI
Rank your gene list by combining layers of evidence in atherosclerosis. This tool allows you to aggregate pre-curated biological data and identify the most credible targets in your list.

# Basic algorithm
Our Basic algorithm has the purpose of scoring each biological layer with user-defined weights.  Internally, the algorithm uses data transformation functions specific to each data type: the raw data is transformed into the same interval of values [0-1] and later aggregated to score every gene independently.

R script funtions: 
