#!/usr/bin/env python3
#
#
# merge HSPs into subject-specific summary statistics
#
# INPUT: output from blastn -fmt 5 (XML)
#
from Bio.Blast import NCBIXML

def parse_blast_output(xml_file):
    print("Open: ", xml_file)
    with open(xml_file, 'r') as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                print(f"# alignment: {alignment.hit_id} len: {alignment.length}" )
                total_length = 0
                identical_matches = 0

                for hsp in alignment.hsps:
                    total_length += hsp.align_length
                    identical_matches += hsp.identities

                identity_percentage = (identical_matches / total_length) * 100 if total_length > 0 else 0

                print(f"Query: {blast_record.query}")
                print(f"Subject: {alignment.hit_def}")
                print(f"Total Length: {total_length}")
                print(f"Identity: {identical_matches}/{total_length} ({identity_percentage:.2f}%)")
                print('-' * 20)

if __name__ == "__main__":
    xml_file = 'blast_output.xml'  # Path to your BLAST XML output
    xml_file = 'results/blastn10_test/a/Keyvirus/JX080302.5.txt'
    parse_blast_output(xml_file)

