from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess
from io import StringIO  # Import StringIO to treat the string as a file-like object

# Set up output file
output_file = open("blast_results.txt", "w")

# Parse human sequences
for record in SeqIO.parse("human.fa", "fasta"):
    human_id = record.id
    sequence = record.seq

    # Run BLAST locally
    result = subprocess.run(
        [
            "/usr/bin/blastn",  # use "blastp" if dealing with protein sequences
            "-query", "-",  # reads sequence from stdin
            "-db", "mouse_db",
            "-outfmt", "5"  # XML format for easy parsing
        ],
        input=str(sequence),
        universal_newlines=True,  # Handle the output as text
        stdout=subprocess.PIPE,   # Captures standard output
        stderr=subprocess.PIPE    # Captures standard error
    )

    # Convert the BLAST output (string) to a file-like object using StringIO
    blast_records = NCBIXML.read(StringIO(result.stdout))

    # Retrieve top hit details
    if blast_records.alignments:
        top_hit = blast_records.alignments[0]
        mouse_id = top_hit.hit_id
        alignment = top_hit.hsps[0].sbjct
        e_value = top_hit.hsps[0].expect
        bitscore = top_hit.hsps[0].bits

        # Write to output file
        output_file.write(f"Human ID: {human_id}, Mouse ID: {mouse_id}\n")
        output_file.write(f"Alignment: {alignment}\n")
        output_file.write(f"E-value: {e_value}, Bitscore: {bitscore}\n\n")
    else:
        output_file.write(f"Human ID: {human_id}, No significant hit found.\n\n")

output_file.close()
