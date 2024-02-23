#!/usr/bin/env python3
#
#
# merge HSPs into subject-specific summary statistics
#
# INPUT: output from blastn -fmt 5 (XML)
#
# Written by ChatGPT v4 
#
# Prompt1: 
#
# please write a short python program using biopython that can read
# the output from command-line blastn and merge all the HSP results
# into a set of per-hit (query, subject pairs), just as the website does. 
# Use existing packages to keep the code concise. 
# 
# Prompt2 
#
# Very good. Instead of just summing the hsp.align_length and
# hsp.identities, do correct overlap calculations to compute what
# percentage of the query sequence is covered by the any hsp. When
# processing and HSP, compute how much of that HSP covers query
# sequence that has not yet been overlapped by an HSP, and only add in
# the percentage of the hsp's identity score that matches the
# percentage of the hsp's length that is not already covered by
# previous hsps.
#
from Bio.Blast import NCBIXML

verbose=False

def parse_blast_output(xml_file):
    print(f"align_summary\tblast_record.id\talignment.hit_id\teffective_length\teffective_identities\tidentity_percentage\tcoverage_percentage")
    with open(xml_file, 'r') as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            query_length = blast_record.query_length
            covered_segments = []  # List to keep track of covered segments of the query

            for alignment in blast_record.alignments:
                print(f"# alignment: {alignment.hit_id} len: {alignment.length}" )
                effective_length = 0
                effective_identities = 0

                for hsp in alignment.hsps:
                    # Calculate non-overlapping part of the current HSP
                    hsp_start, hsp_end = hsp.query_start, hsp.query_end
                    hsp_length = hsp.align_length
                    hsp_identities = hsp.identities

                    # Adjust for overlaps
                    non_overlapping_length = hsp_length
                    for covered_start, covered_end in covered_segments:
                        if hsp_end <= covered_start or hsp_start >= covered_end:
                            continue  # No overlap
                        # Calculate overlap and adjust non_overlapping_length
                        overlap = min(hsp_end, covered_end) - max(hsp_start, covered_start)
                        non_overlapping_length -= overlap

                    # Update effective_length and identities based on non-overlapping portion
                    if non_overlapping_length > 0:
                        percentage_of_hsp_non_overlapped = non_overlapping_length / hsp_length
                        effective_length += non_overlapping_length
                        effective_identities += hsp_identities * percentage_of_hsp_non_overlapped
                        # Add this HSP's range to covered_segments
                        covered_segments.append((hsp_start, hsp_end))

                # Sort covered_segments for future overlap checks
                covered_segments.sort()

                identity_percentage = (effective_identities / effective_length) * 100 if effective_length > 0 else 0
                coverage_percentage = (effective_length / query_length) * 100 if query_length > 0 else 0

                if verbose: 
                    print(f"Query: {blast_record.query}")
                    print(f"Subject: {alignment.hit_def}")
                    print(f"Effective Length (non-overlapping): {effective_length}")
                    print(f"Identity: {effective_identities:.0f}/{effective_length} ({identity_percentage:.2f}%)")
                    print(f"Coverage of Query Sequence: {coverage_percentage:.2f}%")
                    print('-' * 20)

                print(f"align_summary\t{blast_record.query_id}\t{alignment.hit_id}\t{effective_length}\t{effective_identities}\t{identity_percentage}\t{coverage_percentage}")

if __name__ == "__main__":
    xml_file = 'blast_output.xml'  # Path to your BLAST XML output
    xml_file = 'results/blastn10_test/a/Keyvirus/JX080302.5.txt'
    parse_blast_output(xml_file)
