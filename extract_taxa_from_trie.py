
'''
This script takes tries that store NCBI database information on accession number and taxon ids
and iterates through all accessions in a database. For each accession, we get the taxon id
and store a running total of number of accessions and total cumulative sequence length
for each taxon id. This is then stored in a map of the format
{taxid: {'total_length': int, 'num_accessions': int}} and dumped into a pickle file.
Database input should be 'NR' or 'NT'.

Usage: python extract_taxa_from_trie.py {pickle_output_file_name} {database}


'''
import marisa_trie
import re
import pickle
import sys


taxid2accession_trie = marisa_trie.RecordTrie("30p").mmap("marisa_refs/taxid2accession.marisa")
nt_loc_trie = marisa_trie.RecordTrie("QII").mmap('marisa_refs/nt_loc.marisa')
nr_loc_trie = marisa_trie.RecordTrie("QII").mmap('marisa_refs/nr_loc.marisa')
accession_to_taxid_trie = marisa_trie.RecordTrie("Q").mmap('marisa_refs/accession2taxid.marisa')

pickle_filename = sys.argv[1]
db = sys.argv[2]

relevant_trie = nt_loc_trie if db == 'NT' else nr_loc_trie

def populate_taxon_map(taxon_map, this_trie, bad_accessions):
  count = 0
  num_bad_accessions = 0
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

taxon_map = {}
bad_accessions = []
populate_taxon_map(taxon_map, relevant_trie, bad_accessions)
#print(taxon_map)
print(bad_accessions)
with open(pickle_filename, 'wb') as handle:
  pickle.dump(taxon_map, handle, protocol=pickle.HIGHEST_PROTOCOL) #pickle file for 100 taxa is about 2KB


with open('bad_accessions.txt', 'w') as f:
  f.writelines(bad_accessions)
