import pandas as pd
lncrna_df = pd.read_excel("input/Jule_Cedric_Microarray_Reanalysis/reannotation_arraystar/LncRNA Expression Profiling Data.xlsx",skiprows=66)
mrna_df = pd.read_excel("input/Jule_Cedric_Microarray_Reanalysis/reannotation_arraystar/mRNA Expression Profiling Data.xlsx",skiprows=38)
set1 = set(lncrna_df['Transcript_ID'])
set2 = set(lncrna_df['Transcript_ID(hg38)'])

# Overlaps (common elements)
overlaps = set1 & set2  # or: set1.intersection(set2)

# Unique to col1
unique_col1 = set1 - set2

# Unique to col2
unique_col2 = set2 - set1
