import pandas as pd

from workflow.scripts.count import mpranalyze_compiler as compiler


class TestMpranalyzeCompiler:
    def test_get_annot_parses_dna_and_rna_headers(self):
        assert compiler.get_annot("DNA(condition X, replicate 1)") == ("DNA", "X", "1")
        assert compiler.get_annot("RNA(condition Y, replicate 2)") == ("RNA", "Y", "2")

    def test_get_annot_returns_none_for_unmatched_headers(self):
        assert compiler.get_annot("Sequence") == (None, None, None)
        assert compiler.get_annot("label") == (None, None, None)

    def test_generate_annotation_output_repeats_and_numbers_barcodes(self):
        input_frame = pd.DataFrame(
            [
                {"type": "DNA", "condition": "X", "replicate": "1"},
                {"type": "RNA", "condition": "X", "replicate": "1"},
            ]
        )

        output = compiler.generate_annotation_output(input_frame, number_barcodes=2)

        assert list(output["sample"]) == ["DNA_X_1_1", "DNA_X_1_2", "RNA_X_1_1", "RNA_X_1_2"]
        assert list(output["barcode"]) == ["1", "2", "1", "2"]

    def test_generate_count_output_pads_and_flattens_by_barcode(self):
        input_frame = pd.DataFrame(
            [
                {"label": "oligoA", "bc1": 10, "bc2": 20},
                {"label": "oligoA", "bc1": 11, "bc2": 21},
                {"label": "oligoB", "bc1": 12, "bc2": 22},
            ]
        ).set_index("label")

        output = compiler.generate_count_output(input_frame, ["sample_1", "sample_2", "sample_3", "sample_4"], number_barcodes=2)

        assert list(output["seq_id"]) == ["oligoA", "oligoB"]
        assert list(output.loc[output["seq_id"] == "oligoA", ["sample_1", "sample_2", "sample_3", "sample_4"]].iloc[0]) == [10, 11, 20, 21]
        assert list(output.loc[output["seq_id"] == "oligoB", ["sample_1", "sample_2", "sample_3", "sample_4"]].iloc[0]) == [12, 0, 22, 0]

    def test_cli_generates_expected_count_tables(self, cli_runner, minimal_count_input, tmp_path):
        rna_counts_output = tmp_path / "rna_counts.tsv.gz"
        dna_counts_output = tmp_path / "dna_counts.tsv.gz"
        rna_annotation_output = tmp_path / "rna_annot.tsv.gz"
        dna_annotation_output = tmp_path / "dna_annot.tsv.gz"

        result = cli_runner.invoke(
            compiler.cli,
            [
                "--input",
                str(minimal_count_input),
                "--rna-counts-output",
                str(rna_counts_output),
                "--dna-counts-output",
                str(dna_counts_output),
                "--rna-annotation-output",
                str(rna_annotation_output),
                "--dna-annotation-output",
                str(dna_annotation_output),
            ],
        )

        assert result.exit_code == 0, result.output

        rna: pd.DataFrame = pd.read_csv(rna_counts_output, sep="\t").set_index("seq_id")
        dna: pd.DataFrame = pd.read_csv(dna_counts_output, sep="\t").set_index("seq_id")
        cols= ["RNA_X_1_1", "RNA_X_1_2", "RNA_X_2_1", "RNA_X_2_2", "RNA_X_3_1", "RNA_X_3_2"]
        row = rna.loc["oligoA"]
        assert list(row.loc[cols]) == [
            100,
            101,
            200,
            201,
            300,
            301,
        ]
        cols= ["DNA_X_1_1", "DNA_X_1_2", "DNA_X_2_1", "DNA_X_2_2", "DNA_X_3_1", "DNA_X_3_2"]
        row = dna.loc["oligoA"]
        assert list(row.loc[cols]) == [
            10,
            11,
            20,
            21,
            30,
            31,
        ]
