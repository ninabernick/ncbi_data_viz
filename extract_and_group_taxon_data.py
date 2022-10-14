'''
This encapsulates the process of extracting accession data from a trie, fetching taxon data,
cleaning and grouping the data in one script. 
'''
import marisa_trie
import re
import sys
import taxoniq
import csv
import pandas as pd 

# Functions


def populate_taxon_map(this_trie, n=0):
  count = 0
  num_bad_accessions = 0
  bad_accessions = []
  taxon_map = {}
  # keys in the trie are accessions, so iterate over all accessions
  for k in this_trie:
        count +=1
        entry = this_trie[k][0]
        # entry in trie is of the form (seq_offset, header_len, seq_len)
        length = entry[2]
        accession = re.sub("\..*", "", k) #trim . off end
        #get taxon id if possible
        try:  
          taxid = accession_to_taxid_trie[accession][0][0] # this will fail since some records have been removed -- why are they present in nt loc trie?
          # update cumulative length and number of accessions
          if taxid in taxon_map: 
            taxon_map[taxid]["total_length"] = taxon_map[taxid]["total_length"] + length
            taxon_map[taxid]["num_accessions"] = taxon_map[taxid]["num_accessions"] + 1
          else:
            taxon_map[taxid] = {"total_length": length, "num_accessions": 1}
        except:
          num_bad_accessions += 1
          bad_accessions.append(accession + ",")

        if n != 0 and count > n:
          # user set a maximum number of keys to look at
          break
        
  return (taxon_map, bad_accessions)


def create_lineage_array(taxon_map, ranks):
  rows = []
  csv_headers = ranks + ['taxon_id', 'total_length', 'num_accessions']
  rows.append(csv_headers)
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
  return (rows, bad_taxids)

def backfill_data(ranks, df): 
  for rank in ranks:
    df[rank] = df[rank].fillna(f'unknown_{rank}')

  return df

# Script

store_bad_data = False

# fetch necessary tries
taxid2accession_trie = marisa_trie.RecordTrie("30p").mmap("marisa_refs/taxid2accession.marisa")
nt_loc_trie = marisa_trie.RecordTrie("QII").mmap('marisa_refs/nt_loc.marisa')
nr_loc_trie = marisa_trie.RecordTrie("QII").mmap('marisa_refs/nr_loc.marisa')
accession_to_taxid_trie = marisa_trie.RecordTrie("Q").mmap('marisa_refs/accession2taxid.marisa')

db = sys.argv[1]
csv_filename = sys.argv[2]
output_csv_name = sys.argv[3]

relevant_trie = nt_loc_trie if db == 'NT' else nr_loc_trie


# populate a taxon map of the form {taxid: {'total_length': int, 'num_accessions': int}}
(taxon_map, bad_accessions) = populate_taxon_map(relevant_trie)

# store accessions that didn't map to taxon ids
if store_bad_data:
  with open('invalid_accessions.txt', 'w') as f:
    f.writelines(bad_accessions)
  
# create 2D array with taxon lineage info, first row is column headers
ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

data, bad_taxids = create_lineage_array(taxon_map, ranks)

# save intermediary csv before continuing
with open(csv_filename, 'w', newline='') as f:
  writer = csv.writer(f, delimiter=",")
  writer.writerows(data) 

if store_bad_data:
  with open("bad_taxa.txt", 'w', newline='') as f:
    f.writelines(bad_taxids)


# clean and group data
all_data = pd.read_csv(csv_filename)
# group taxa under "other sequences" into one row
other_sequences = all_data[all_data.superkingdom.isna()].sort_values(by='num_accessions', ascending=False)
length_sum = other_sequences.total_length.sum()
accession_sum = other_sequences.num_accessions.sum()
all_data.dropna(subset='superkingdom', inplace=True)
data = ["other_sequences"] * len(ranks) + [-1, length_sum, accession_sum]
df2 = pd.DataFrame([data], columns = all_data.columns)
all_data = pd.concat([df2, all_data], ignore_index=True)

# Bacteria and Archaea are missing kingdoms. Fill in the kingdom level with the superkingdom.
all_data.loc[all_data.superkingdom == 'Archaea', 'kingdom']= all_data[all_data.superkingdom == 'Archaea'].kingdom.fillna(all_data.superkingdom)
all_data.loc[all_data.superkingdom == 'Bacteria', 'kingdom'] = all_data[all_data.superkingdom == 'Bacteria'].kingdom.fillna(all_data.superkingdom)

# fill missing ranks to get rid of NAs
# all data has superkingdom at least
all_data = backfill_data(ranks, all_data)

groupby_col = 'order'
df_grouped = all_data.groupby(ranks[0: ranks.index(groupby_col)+1], as_index=False).agg(
  {
    'family': 'unique',
    'genus': 'unique',
    'species': 'unique',
    'taxon_id': 'unique',
    'total_length': 'sum', 
    'num_accessions': 'sum'
  })

df_grouped.to_csv(output_csv_name)