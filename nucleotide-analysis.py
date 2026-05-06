from Bio import Entrez, SeqIO
import pandas as pd


# Your NCBI registered email — required by Entrez
Entrez.email = "ketz.bee@gmail.com"

# Accession IDs to analyse — edit this list to change your targets
ACCESSION_IDS = ["NM_002232.5", "NM_000546.6", "NM_007294.3"]


def fetch_genbank_records(accession_ids):
    """
    Fetch GenBank records from NCBI for a list of accession IDs.

    Args:
        accession_ids (list): List of NCBI accession ID strings.

    Returns:
        list: List of SeqRecord objects in GenBank format.
    """
    handle = Entrez.efetch(
        db="nucleotide",
        id=",".join(accession_ids),
        rettype="gb",
        retmode="text"
    )
    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    return records


def compute_nucleotide_stats(record):
    """
    Extract CDS from a GenBank record and compute nucleotide counts and ratios.

    Args:
        record: A Biopython SeqRecord object in GenBank format.

    Returns:
        dict of nucleotide counts and AT/GC ratios, or None if CDS is invalid.
    """
    for feature in record.features:
        if feature.type == "CDS":
            cds_sequence = feature.extract(record.seq)

            num_a = cds_sequence.count("A")
            num_t = cds_sequence.count("T")
            num_g = cds_sequence.count("G")
            num_c = cds_sequence.count("C")

            # Skip if any nucleotide is absent — likely an invalid sequence
            if num_a == 0 or num_t == 0 or num_g == 0 or num_c == 0:
                print(f"Skipping {record.name} — invalid or incomplete sequence")
                return None

            return {
                "A": num_a,
                "T": num_t,
                "G": num_g,
                "C": num_c,
                "AT_ratio": round(num_a / num_t, 3),
                "GC_ratio": round(num_g / num_c, 3)
            }


def build_results_table(records):
    """
    Process a list of GenBank records and return a pandas DataFrame
    with nucleotide stats for each accession ID.

    Args:
        records (list): List of Biopython SeqRecord objects.

    Returns:
        pd.DataFrame: Table with one row per accession ID.
    """
    data = {}

    for record in records:
        stats = compute_nucleotide_stats(record)
        if stats is not None:
            data[record.name] = stats

    df = pd.DataFrame(data).T
    df.index.name = "Accession ID"
    return df


def main():
    print("Fetching records from NCBI...")
    records = fetch_genbank_records(ACCESSION_IDS)

    print("Computing nucleotide statistics...\n")
    results = build_results_table(records)

    print(results)


main()
