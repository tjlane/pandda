import numpy as np

from giant.structure.sequence import (
    align_sequences_default,
    pairwise_sequence_identity,
    )

seq_a = 'SMSVKKPKRDDSKDLALCSMILTEMETHEDAWPFLLPVNLKLVPGYKKVIKKPMDFSTIREKLSSGQYPNLETFALDVRLVFDNCETFNEDDSDIGRAGHNMRKYFEKKWTDTFK'
seq_b = 'SMSVKKPKRDDSKDLALCSMILTEMETHEDAWPFLLPVNLKLVPGYKKVIKKPMDFSTIREKLSSGQYPNLETFALDVRLVFDNCETFNEDDSDIGRAGHNMRKYFEKKW' # Truncated sequence
seq_c = seq_a.replace('VNL','AAA') # Mutated sequence
seq_d = 'SMSVKKPKRDDSKDLALCSMILTEMETHEDAWPFLLPVAQNPNCNIMIFHPTKEEFNDFDKYIAYMESQGAHRAGLAKIIPPKEWKARETYDNISEILIATPLQQVASGRAGVFT' # Spliced with unrelated (half overlaps)
seq_e = 'AQNPNCNIMIFHPTKEEFNDFDKYIAYMESQGAHRAGLAKIIPPKEWKARETYDNISEILIATPLQQVASGRAGVFTQYHKKKKAMTVGEYRHLANSKKYQTPPHQNFEDLERKYWKNRIYNSPIYGADISGSLFDENTKQWNLG' # Unrelated sequence

set_1 = [seq_a, seq_c, seq_d, seq_e]
set_2 = [seq_b, seq_c, seq_d]


def test_align_sequences_default_identical():

    ali = align_sequences_default(seq_a, seq_a)

    assert ali.a == seq_a
    assert ali.b == seq_a
    assert ali.calculate_sequence_identity() == 1.0
    assert ali.matches() == '|'*len(seq_a)

def test_align_sequences_default_truncated():

    ali = align_sequences_default(seq_a, seq_b)

    assert ali.a == seq_a
    assert ali.b == seq_b+'-'*5
    assert ali.calculate_sequence_identity() == 110/115.
    assert ali.matches() == '|'*len(seq_b)+' '*5

def test_align_sequences_default_mutation():

    ali = align_sequences_default(seq_a, seq_c)

    assert ali.a == seq_a
    assert ali.b == seq_c
    assert ali.calculate_sequence_identity() == 112/115.
    assert ali.matches() == '|'*37+' '*3+'|'*75

def test_pairwise_sequence_identity():

    arr = pairwise_sequence_identity(
        seqs_1 = set_1,
        seqs_2 = set_2,
        min_alignment = 0.0,
        seq_identity_threshold = None,
        )

    assert np.array_equal(
        arr.round(5),
        np.array([
            [ 110/115.   , 112/115.   , 42/115.   ],
            [ 107/115.   , 1.         , 42/115.   ],
            [ 42/115.    , 42/115.    , 1.        ],
            [ 19/145.    , 15/145.    , 77/145.   ],
            ]).round(5),
        )

def test_pairwise_sequence_identity_frac_sequence_alignment():

    arr = pairwise_sequence_identity(
        seqs_1 = set_1,
        seqs_2 = set_2,
        min_alignment = 0.60,
        seq_identity_threshold = None,
        )

    assert np.array_equal(
        arr.round(5),
        np.array([
            [ 110/115.   , 112/115.   , 0.        ],
            [ 107/115.   , 1.         , 0.        ],
            [ 0.         , 0.         , 1.        ],
            [ 0.         , 0.         , 77/145.   ],
            ]).round(5),
        )

def test_pairwise_sequence_identity_num_sequence_alignment():

    arr = pairwise_sequence_identity(
        seqs_1 = set_1,
        seqs_2 = set_2,
        min_alignment = 60,
        seq_identity_threshold = None,
        )

    assert np.array_equal(
        arr.round(5),
        np.array([
            [ 110/115.   , 112/115.   , 0.        ],
            [ 107/115.   , 1.         , 0.        ],
            [ 0.         , 0.         , 1.        ],
            [ 0.         , 0.         , 77/145.   ],
            ]).round(5),
        )

def test_pairwise_sequence_identity_seq_identity_threshold():

    arr = pairwise_sequence_identity(
        seqs_1 = set_1,
        seqs_2 = set_2,
        min_alignment = 30,
        seq_identity_threshold = 0.90,
        )

    assert np.array_equal(
        arr,
        np.array([
            [ 1 , 1 , 0 ],
            [ 1 , 1 , 0 ],
            [ 0 , 0 , 1 ],
            [ 0 , 0 , 0 ],
            ]),
        )
