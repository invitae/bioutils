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
        ("ATTAATCCC", TranslationTable.vertebrate_mitochondrial, True, "MNP"),
        ("ATTATTAATCCC", TranslationTable.vertebrate_mitochondrial, True, "MINP"),
        ("ATTAATCCC", TranslationTable.vertebrate_mitochondrial, False, "INP"),
        ("ATTAATCCC", TranslationTable.standard, True, "INP"),
        ("ATTAATCCC", TranslationTable.standard, False, "INP"),
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
    assert (translate_cds(sequence, full_codons=False) == translated_sequence)


def test_translate_cds_full_codons_true():
    assert (translate_cds("TTT", full_codons=True) == "F")

    with pytest.raises(ValueError):
        translate_cds("TT", full_codons=True)
        translate_cds("TTAAA", full_codons=True)