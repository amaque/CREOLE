"""
3A test -> AAA -> Arange - Act - Assets
"""
from pytest import raises

from creole.creole_prot import orf_prediction as orf_p


def test_nucleotide_sequences_must_be_capital_letters():
    # Arange
    id_seq: str = '>read_001'
    sequence: str = 'atgctgcttcgtatagctagtactagct'
    min_orf_len: int = 75

    # Act
    result = orf_p.OrfsPrediction(id_seq, sequence, min_orf_len)
    # Assert
    assert result


def test_raise_an_error_if_nucleotide_is_not_valid():
    id_seq: str = '>read_001'
    sequence: str = 'atgctRcttcgtaVagctagtactagct'
    min_orf_len: int = 75
    seq_ype = ['A', 'C', 'G', 'T']

    error_msg = (
        f'Provided data does not seem to be a correct {seq_ype} sequence'
    )
    # with raises() as error:
    orf_p.OrfsPrediction(id_seq, sequence, min_orf_len)

    # assert error_msg == error.value.args[0]
