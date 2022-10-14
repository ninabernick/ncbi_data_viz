'''
This file reads in a taxon map of the form 
{taxid: {'total_length': int, 'num_accessions': int}}
then tries to use taxoniq to find lineage information for each taxon.

Results are output into a csv with headers:
'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'taxon_id', 'total_length', 'num_accessions'

Missing values in lineage information are filled with empty strings. Processing of these empty values is done in the cleaning
and grouping step of the next file. 

Usage: python write_csv.py {taxon_map_pickle_filename} {output_csv_filename}
'''
import sys
import pickle
import taxoniq
import csv

pickle_filename = sys.argv[1]
csv_filename = sys.argv[2]

with open(pickle_filename, 'rb') as f:
  taxon_map = pickle.load(f)

rows = []
ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
csv_headers = ranks + ['taxon_id', 'total_length', 'num_accessions']
bad_taxids = []
for taxid, data in taxon_map.items():
  taxon_row = []
  lineage = {}

  try:
    taxon = taxoniq.Taxon(taxid)
    for t in taxon.ranked_lineage:
      lineage[t.rank.name] = t.scientific_name
  except: 
    bad_taxids.append(str(taxid) + ',')

  for r in ranks:
    lineage_member = lineage.get(r, "")
    taxon_row.append(lineage_member)

  taxon_row.append(taxid)
  taxon_row.append(data.get("total_length", ""))
  taxon_row.append(data.get("num_accessions", ""))

  rows.append(taxon_row)

print("bad taxon ids")
print(bad_taxids)
with open(csv_filename, 'w', newline='') as f:
  writer = csv.writer(f, delimiter=",")
  writer.writerow(csv_headers)
  writer.writerows(rows) 

with open("bad_taxa.txt", 'w', newline='') as f:
  f.writelines(bad_taxids)


