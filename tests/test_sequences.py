import pytest

from bioutils.sequences import TranslationTable, translate_cds


def test_translate_examples():
    """test for standard translation table"""

    assert translate_cds("ATGCGA") == "MR"
    assert translate_cds("AUGCGA") == "MR"
    assert translate_cds(None) is None
    assert translate_cds("") == ""
    with pytest.raises(ValueError):
        translate_cds("AUGCG")

    assert translate_cds("AUGCG", full_codons=False) == "M*"
    assert translate_cds("ATGTAN") == "MX"
    assert translate_cds("CCN") == "P"
    assert translate_cds("TRA") == "*"
    assert translate_cds("TTNTA", full_codons=False) == "X*"
    assert translate_cds("CTB") == "L"
    assert translate_cds("AGM") == "X"
    assert translate_cds("GAS") == "X"
    assert translate_cds("CUN") == "L"
    with pytest.raises(ValueError):
        translate_cds("AUGCGQ")


def test_translate_selenoproteins():
    """unit test for sec codon"""
    assert translate_cds("AUGTGATAA") == "M**"
    assert translate_cds("AUGTGATAA", translation_table=TranslationTable.standard) == "M**"
    assert translate_cds("AUGTGATAA", translation_table=TranslationTable.selenocysteine) == "MU*"
    assert (
        translate_cds(
            "AUGTGATA",
            translation_table=TranslationTable.selenocysteine,
            full_codons=False,
        )
        == "MU*"
    )

    with pytest.raises(ValueError):
        translate_cds("AUGTGATA", translation_table=TranslationTable.selenocysteine)


def test_translate_vertebrate_mitochondrial():
    """unit test for vertebrate mitochondrial codons"""
    assert translate_cds("AUGTGATAA") == "M**"
    assert translate_cds("ATATGAAGGAGA", translation_table=TranslationTable.vertebrate_mitochondrial) == "MW**"
    assert (
        translate_cds(
            "ATAAG",
            translation_table=TranslationTable.vertebrate_mitochondrial,
            full_codons=False,
        )
        == "M*"
    )

    with pytest.raises(ValueError):
        translate_cds("ATAAG", translation_table=TranslationTable.vertebrate_mitochondrial)


@pytest.mark.parametrize(
    "sequence, exception_map, translated_sequence",
    (
        ("ATGATGATG", {3: "U"}, "MUM"),
        ("ATGATGATG", {3: "U", 6: "U"}, "MUU"),
        ("ATGATGATG", {6: "*"}, "MM*"),
        ("ATGATGAT", {6: "*"}, "MM*"),
        ("ATGACTATG", {}, "MTM"),
        ("ATGACTATG", None, "MTM"),
    ),
)
def test_translate_cds_w_exceptions(sequence, exception_map, translated_sequence):
    assert translate_cds(sequence, full_codons=False, ter_symbol="", exception_map=exception_map) == translated_sequence


@pytest.mark.parametrize(
    "sequence, translation_table, starts_at_first_codon, translated_sequence",
    (
        pytest.param(
            "ATTAATCCC",
            TranslationTable.vertebrate_mitochondrial,
            True,
            "MNP",
            id="mitochondrial alternate start codon ATT when at the first position",
        ),
        pytest.param(
            "ATTATTAATCCC",
            TranslationTable.vertebrate_mitochondrial,
            True,
            "MINP",
            id="mitochondrial alternate start ATT only used only at first position",
        ),
        pytest.param(
            "ATTAATCCC",
            TranslationTable.vertebrate_mitochondrial,
            False,
            "INP",
            id="mitochondrial alternate start ATT not used if first position is False",
        ),
        pytest.param(
            "ATTAATCCC",
            TranslationTable.standard,
            True,
            "INP",
            id="ATT not an alternate start for human nuclear codons",
        ),
        pytest.param(
            "ATTAATCCC",
            TranslationTable.standard,
            False,
            "INP",
            id="start position is False so alternate starts are not used",
        ),
        pytest.param(
            "ATGAATCCC",
            TranslationTable.standard,
            True,
            "MNP",
            id="ATG is a standard start codon for nuclear codons",
        ),
        pytest.param(
            "ATGAATCCC",
            TranslationTable.standard,
            False,
            "MNP",
            id="ATG is a standard start codon for nuclear even if not first position",
        ),
        pytest.param(
            "CTGAATCCC",
            TranslationTable.standard,
            True,
            "MNP",
            id="CTG is a valid alternate start codon for nuclear codons if at first position",
        ),
        pytest.param(
            "CTGAATCCC",
            TranslationTable.vertebrate_mitochondrial,
            True,
            "LNP",
            id="CTG is not a valid alternate start codon for mitochondrial codons",
        ),
        pytest.param(
            "CTGAATCCC",
            TranslationTable.standard,
            False,
            "LNP",
            id="CTG is a L codon for nuclear codons if not at first position",
        ),
        pytest.param(
            "TTGAATCCC",
            TranslationTable.standard,
            True,
            "MNP",
            id="TTG is a valid alternate start codon for nuclear codons if at first position",
        ),
        pytest.param(
            "TTGAATCCC",
            TranslationTable.vertebrate_mitochondrial,
            True,
            "LNP",
            id="TTG is not a valid alternate start codon for mitochondral codons",
        ),
        pytest.param(
            "TTGAATCCC",
            TranslationTable.standard,
            False,
            "LNP",
            id="TTG is a L codon for nuclear codons if not at first position",
        ),
        pytest.param(
            "GTGAATCCC",
            TranslationTable.standard,
            True,
            "VNP",
            id="GTG is not a valid alternate start codon for nuclear codons",
        ),
    ),
)
def test_translate_cds_mito_alts(sequence, translation_table, starts_at_first_codon, translated_sequence):
    assert (
        translate_cds(
            sequence,
            full_codons=False,
            ter_symbol="",
            translation_table=translation_table,
            starts_at_first_codon=starts_at_first_codon,
        )
        == translated_sequence
    )


@pytest.mark.parametrize(
    "sequence, translated_sequence",
    (
        ("ATTATTA", "II*"),
        ("ATTATT", "II"),
        ("ATTAT", "I*"),
        ("ATTA", "I*"),
        ("GG", "*"),
        ("G", "*"),
    ),
)
def test_translate_cds_full_codons_false(sequence, translated_sequence):
    assert translate_cds(sequence, full_codons=False) == translated_sequence


def test_translate_cds_full_codons_true():
    assert translate_cds("TTT", full_codons=True) == "F"

    with pytest.raises(ValueError):
        translate_cds("TT", full_codons=True)
