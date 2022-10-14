# NCBI Data Visualization
Project to visualize NCBI database contents as part of CZI Science Make-a-Thon. CZID maps sequences to NCBI database sequences, so of particular interest is determining total sequence length for each taxon or higher taxonomy rank.

The goal of these scripts is to extract NCBI data from the stored trie format into a CSV that provides taxonomy lineage for a given taxon and aggregates the total length of sequences in the NCBI database and total number of acccessions.
## Running Scripts:

To generate data in one step, use `extract_and_group_taxon_data.py` as follows:

`>> python extract_and_group_taxon_data.py {NT | NR} {intermediary_csv_filename} {output_csv_filename}`

where NT or NR is the desired database, and the csv filenames store intermediary output and final output (taxon information before and after cleaning and grouping, respectively).

To generate data step by step, use the other three scripts:

```
>> python extract_taxa_from_trie.py {taxon_map_pickle_filename} {database}
>> python write_csv.py {taxon_map_pickle_filename} {output_csv_filename}
```
and then run the .ipynb file, setting desired level of grouping and appropriate file names. 
