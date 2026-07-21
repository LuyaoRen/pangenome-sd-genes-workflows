import os
import argparse
from Bio import SeqIO
import subprocess
from io import StringIO

def map_sequences_to_single_paf(shared_fasta, target_file, combined_output_file):
    """Map each sequence in the multi-FASTA input file to the target using minimap2, appending results to a single output file."""
    
    # Open the combined output file in write mode
    with open(combined_output_file, "w") as out_f:
        for record in SeqIO.parse(target_file, "fasta"):
            # Prepare a FASTA formatted string for the sequence
            fasta_string = f">{record.id}\n{record.seq}\n"
            
            # Use StringIO to simulate a file object for the sequence
            with StringIO(fasta_string) as fasta_io:
                # Define minimap2 command with options as separate items
                minimap2_command = [
                    "minimap2", "-x", "asm10", "-c", "--secondary=yes", 
                    "-p", "0.3", "-N", "10000", "--eqx", "-t", "4", 
                    "-r", "500", "-K", "500M", "-" , shared_fasta
                ]
                
                # Run minimap2, providing the sequence via stdin
                result = subprocess.run(minimap2_command, input=fasta_io.read(), text=True, stdout=subprocess.PIPE)
                
                # Write the result to the combined output file
                out_f.write(result.stdout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Get shared bed region for each sequence in a FASTA file with multiple sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-s", "--shared_fa", type=str, help="shared fasta file", required=True)
    parser.add_argument("-t", "--target_fa", type=str, help="target fasta file", required=True)
    parser.add_argument("-p", "--paf_file", type=str, help="output PAF file", required=True)

    args = parser.parse_args()

    # Set paths from arguments
    shared_fasta = args.shared_fa
    target_file = args.target_fa  
    combined_output_file = args.paf_file

    # Run the mapping function
    map_sequences_to_single_paf(shared_fasta, target_file, combined_output_file)
