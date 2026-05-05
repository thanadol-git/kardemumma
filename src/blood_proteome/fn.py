import requests
import sys
import json
import os
import time

import pridepy

WEBSITE_API = "https://rest.uniprot.org/"
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"


def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)
    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()
    return response


def get_protein_sequences_batch(protein_ids, batch_size=50):
    """Retrieve protein sequences from UniProt in batches to avoid URL length limits."""
    all_sequences = []
    clean_protein_ids = [pid.strip() for pid in protein_ids if pid and pid.strip()]

    print(f"Total proteins to query: {len(clean_protein_ids)}")

    for i in range(0, len(clean_protein_ids), batch_size):
        batch = clean_protein_ids[i:i + batch_size]
        n_batches = (len(clean_protein_ids) - 1) // batch_size + 1
        print(f"Processing batch {i // batch_size + 1}/{n_batches} ({len(batch)} proteins)")

        accessions = " OR ".join([f"accession:{acc}" for acc in batch])
        fields = "accession,sequence,protein_name,gene_names,organism_name,organism_id,reviewed,id"
        query_url = f"{WEBSITE_API}uniprotkb/search?query=({accessions})&fields={fields}&format=json"

        try:
            response = get_url(query_url)
            data = response.json()
            if 'results' in data:
                all_sequences.extend(data['results'])
                print(f"  Retrieved {len(data['results'])} sequences")
            else:
                print(f"  No results in batch {i // batch_size + 1}")
        except Exception as e:
            print(f"  Error in batch {i // batch_size + 1}: {e}")
            continue

        time.sleep(0.1)

    print(f"Total sequences retrieved: {len(all_sequences)}")
    return all_sequences


def get_swissprot_sequences_batch(protein_ids, batch_size=50):
    """Retrieve only Swiss-Prot reviewed protein sequences from UniProt in batches."""
    all_sequences = []
    clean_protein_ids = [pid.strip() for pid in protein_ids if pid and pid.strip()]

    print(f"Total proteins to query: {len(clean_protein_ids)}")
    print("Filtering for Swiss-Prot reviewed entries only...")

    for i in range(0, len(clean_protein_ids), batch_size):
        batch = clean_protein_ids[i:i + batch_size]
        n_batches = (len(clean_protein_ids) - 1) // batch_size + 1
        print(f"Processing batch {i // batch_size + 1}/{n_batches} ({len(batch)} proteins)")

        accessions = " OR ".join([f"accession:{acc}" for acc in batch])
        full_query = f"({accessions}) AND reviewed:true"
        fields = "accession,sequence,protein_name,gene_names,organism_name,organism_id,reviewed,id"
        query_url = f"{WEBSITE_API}uniprotkb/search?query={full_query}&fields={fields}&format=json"

        try:
            response = get_url(query_url)
            data = response.json()
            if 'results' in data:
                reviewed_entries = [
                    e for e in data['results'] if 'Swiss-Prot' in e.get('entryType', '')
                ]
                all_sequences.extend(reviewed_entries)
                print(f"  Retrieved {len(reviewed_entries)} Swiss-Prot sequences "
                      f"(out of {len(data['results'])} total)")
            else:
                print(f"  No results in batch {i // batch_size + 1}")
        except Exception as e:
            print(f"  Error in batch {i // batch_size + 1}: {e}")
            continue

        time.sleep(0.1)

    print(f"Total Swiss-Prot sequences retrieved: {len(all_sequences)}")
    return all_sequences


def query_human_proteome(reviewed_only=True, batch_size=500):
    """
    Query the entire human proteome from UniProt.

    Parameters
    ----------
    reviewed_only : bool
        If True, return only Swiss-Prot reviewed proteins (default: True).
    batch_size : int
        Proteins per batch (default: 500, unused for stream endpoint).

    Returns
    -------
    list
        Protein entries with sequences and metadata.
    """
    print("=== Querying Human Proteome from UniProt ===")

    if reviewed_only:
        query = "organism_id:9606 AND reviewed:true"
        print("Retrieving Swiss-Prot reviewed human proteins only...")
    else:
        query = "organism_id:9606"
        print("Retrieving all human proteins (Swiss-Prot + TrEMBL)...")

    count_url = f"{WEBSITE_API}uniprotkb/search?query={query}&size=0"
    try:
        response = get_url(count_url)
        total_count = int(response.headers.get('x-total-results', 0))
        print(f"Total human proteins available: {total_count:,}")
    except Exception as e:
        print(f"Error getting count: {e}")
        return []

    if total_count == 0:
        print("No proteins found!")
        return []

    fields = "accession,sequence,protein_name,gene_names,organism_name,organism_id,reviewed,id"
    stream_url = f"{WEBSITE_API}uniprotkb/stream?query={query}&fields={fields}&format=json"
    print(f"Starting download of {total_count:,} proteins...")

    try:
        response = get_url(stream_url)
        data = response.json()
        all_sequences = data.get('results', [])
        print(f"Successfully retrieved {len(all_sequences):,} human protein sequences")

        if reviewed_only:
            reviewed_count = len(all_sequences)
            unreviewed_count = 0
        else:
            reviewed_count = sum(1 for s in all_sequences if 'Swiss-Prot' in s.get('entryType', ''))
            unreviewed_count = len(all_sequences) - reviewed_count

        print(f"Reviewed (Swiss-Prot): {reviewed_count:,}")
        print(f"Unreviewed (TrEMBL): {unreviewed_count:,}")

        if all_sequences:
            first = all_sequences[0]
            print(f"\nExample entry:")
            print(f"  Accession: {first.get('primaryAccession', 'N/A')}")
            print(f"  Entry type: {first.get('entryType', 'N/A')}")
            print(f"  Sequence length: {len(first.get('sequence', {}).get('value', ''))}")

        return all_sequences

    except Exception as e:
        print(f"Error retrieving human proteome: {e}")
        return []


def create_human_proteome_fasta(reviewed_only=True, file_loc=''):
    """
    Download the human proteome from UniProt and write it to a FASTA file.

    Parameters
    ----------
    reviewed_only : bool
        If True, include only Swiss-Prot reviewed proteins (default: True).
    file_loc : str
        Directory to save the FASTA file (default: current directory).

    Returns
    -------
    str or None
        Path to the created FASTA file, or None on failure.
    """
    file_loc = os.path.normpath(str(file_loc).strip()) if str(file_loc).strip() else ''
    if file_loc:
        try:
            os.makedirs(file_loc, exist_ok=True)
        except Exception as e:
            print(f"Warning: could not create directory '{file_loc}': {e}. Writing to current directory.")
            file_loc = ''

    print("=== Creating Human Proteome FASTA File ===")
    human_proteins = query_human_proteome(reviewed_only=reviewed_only)

    if not human_proteins:
        print("No proteins retrieved. Aborting FASTA creation.")
        return None

    today = time.strftime("%y%m%d")
    base_name = (
        f"Human_proteome_SwissProt_{today}.fasta" if reviewed_only
        else f"Human_proteome_All_{today}.fasta"
    )
    fasta_filename = os.path.join(file_loc, base_name) if file_loc else base_name
    print(f"Creating FASTA file: {fasta_filename}")

    try:
        with open(fasta_filename, "w") as fasta_file:
            for i, entry in enumerate(human_proteins):
                if i % 1000 == 0:
                    print(f"  Processing protein {i + 1:,}/{len(human_proteins):,}")

                prot_name = "Unknown protein"
                desc = entry.get('proteinDescription', {})
                if 'recommendedName' in desc and desc['recommendedName']:
                    full_name = desc['recommendedName'].get('fullName', '')
                    prot_name = full_name.get('value', prot_name) if isinstance(full_name, dict) else full_name or prot_name
                elif 'submittedName' in desc and desc['submittedName']:
                    submitted = desc['submittedName']
                    if isinstance(submitted, list) and submitted:
                        full_name = submitted[0].get('fullName', '')
                        prot_name = full_name.get('value', prot_name) if isinstance(full_name, dict) else full_name or prot_name

                accession = entry.get('primaryAccession', 'UNKNOWN')
                sequence = entry.get('sequence', {}).get('value', '')
                if not sequence:
                    continue

                gene_name = ''
                genes = entry.get('genes', [])
                if isinstance(genes, list) and genes:
                    gene_info = genes[0]
                    if isinstance(gene_info, dict):
                        gene_name = gene_info.get('geneName', {}).get('value', '') or gene_info.get('geneName', '')

                reviewed = entry.get('reviewed', None)
                is_reviewed = reviewed if isinstance(reviewed, bool) else 'Swiss-Prot' in entry.get('entryType', '')
                db_type = "sp" if is_reviewed else "tr"
                entry_name = str(entry.get('uniProtkbId') or entry.get('uniProtKBId') or f"{accession}_HUMAN")
                prot_name = prot_name.replace('\n', ' ').replace('\r', ' ').strip()

                taxid = entry.get('organism', {}).get('taxonId') if isinstance(entry.get('organism'), dict) else None
                pe = '1' if is_reviewed else '2'
                sv = str(entry.get('version', 1))

                header_parts = [f"{db_type}|{accession}|{entry_name}", prot_name, "OS=Homo sapiens",
                                f"OX={taxid or 9606}"]
                if gene_name:
                    header_parts.append(f"GN={gene_name}")
                header_parts += [f"PE={pe}", f"SV={sv}"]

                fasta_file.write(f">{' '.join(header_parts)}\n")
                for j in range(0, len(sequence), 60):
                    fasta_file.write(sequence[j:j + 60] + "\n")

        file_size = os.path.getsize(fasta_filename)
        print(f"\nFASTA file created: {fasta_filename}")
        print(f"Proteins: {len(human_proteins):,}  |  Size: {file_size / (1024 * 1024):.1f} MB")
        return fasta_filename

    except Exception as e:
        print(f"Error creating FASTA file: {e}")
        return None


def validate_fasta(filename):
    """Print summary statistics and the first five headers of a FASTA file."""
    try:
        protein_count = 0
        sequence_lengths = []
        current_sequence = ""

        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_sequence:
                        sequence_lengths.append(len(current_sequence))
                        current_sequence = ""
                    protein_count += 1
                else:
                    current_sequence += line

        if current_sequence:
            sequence_lengths.append(len(current_sequence))

        print(f"FASTA validation: {filename}")
        print(f"  Proteins:            {protein_count}")
        print(f"  Avg sequence length: {sum(sequence_lengths) / len(sequence_lengths):.1f} aa")
        print(f"  Shortest:            {min(sequence_lengths)} aa")
        print(f"  Longest:             {max(sequence_lengths)} aa")
        print(f"  Total amino acids:   {sum(sequence_lengths):,}")

        print("\nFirst 5 headers:")
        header_count = 0
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('>') and header_count < 5:
                    print(f"  {line.strip()}")
                    header_count += 1
                if header_count == 5:
                    break

    except Exception as e:
        print(f"Error validating FASTA file: {e}")
