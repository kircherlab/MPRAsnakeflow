import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts", "count"))

def main():
    import tempfile

    import pandas as pd

    import mpranalyze_compiler as c

    here = os.path.dirname(__file__)
    with tempfile.TemporaryDirectory() as tmp:
        rna_counts_output = os.path.join(tmp, "rna_counts.tsv.gz")
        dna_counts_output = os.path.join(tmp, "dna_counts.tsv.gz")
        rna_annot_output = os.path.join(tmp, "rna_annot.tsv.gz")
        dna_annot_output = os.path.join(tmp, "dna_annot.tsv.gz")

        c.cli.callback(
            input_file=os.path.join(here, "minimal_test_input.tsv"),
            rna_counts_output_file=rna_counts_output,
            dna_counts_output_file=dna_counts_output,
            rna_annotation_output_file=rna_annot_output,
            dna_annotation_output_file=dna_annot_output,
        )

        rna = pd.read_csv(rna_counts_output, sep="\t").set_index("seq_id")
        dna = pd.read_csv(dna_counts_output, sep="\t").set_index("seq_id")

        assert list(rna.loc["oligoA", ["RNA_X_1_1", "RNA_X_1_2", "RNA_X_2_1", "RNA_X_2_2", "RNA_X_3_1", "RNA_X_3_2"]]) == [100, 101, 200, 201, 300, 301]
        assert list(dna.loc["oligoA", ["DNA_X_1_1", "DNA_X_1_2", "DNA_X_2_1", "DNA_X_2_2", "DNA_X_3_1", "DNA_X_3_2"]]) == [10, 11, 20, 21, 30, 31]

if __name__=="__main__":
    main()
