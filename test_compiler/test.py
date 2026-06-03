import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts", "count"))

def main():
    import mpranalyze_compiler as c
    here = os.path.dirname(__file__)
    c.cli.callback(
        input_file=os.path.join(here, "minimal_test_input.tsv"),
        rna_counts_output_file=os.path.join(here, "rna_counts.tsv.gz"),
        dna_counts_output_file=os.path.join(here, "dna_counts.tsv.gz"),
        rna_annotation_output_file=os.path.join(here, "rna_annot.tsv.gz"),
        dna_annotation_output_file=os.path.join(here, "dna_annot.tsv.gz"),
    )

if __name__=="__main__":
    main()
