import pytest
from normalizer.normalizer import mutalyzer3
from pathlib import Path


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), 'r') as file:
        content = file.read()
    return content


def fetch_annotation(reference_id, reference_type=None):
    return _get_content('data/' + reference_id + '.gff3'), 'gff', 'ncbi'


def fetch_sequence(reference_id, reference_source=None):
    return _get_content('data/' + reference_id + '.sequence')


@pytest.mark.parametrize(
    'hgvs_description, normalized_description',
    [('NG_012337.1:g.4delins7_31',
      'NG_012337.1:g.4delins7_31'),
     ('NG_012337.1:g.26_31del',
      'NG_012337.1:g.29_34del'),
     ('NG_012337.1:g.4delins7_50',
      'NG_012337.1:g.[3_4insGGTT;4_5ins12_50]'),
     ('NG_012337.1:g.26_31del',
      'NG_012337.1:g.29_34del'),
     ('NG_012337.1:g.4A>T',
      'NG_012337.1:g.4C>T'),
     ('NG_012337.1:g.100_200delins100_101',
      'NG_012337.1:g.102_200del'),
     ('NG_017013.2:g.19258dup',
      'NG_017013.2:g.19258dup'),
     ('NG_017013.2:g.17471_17471del',
      'NG_017013.2:g.17471del'),
     ('NG_017013.2:g.18748_18750delinsCAT',
      'NG_017013.2:g.18749G>A'),
     ('NG_017013.2:g.17394_17395insC',
      'NG_017013.2:g.17394dup'),
     ('NG_012337.1(NM_003002.2):c.274G>T',
      'NG_012337.1:g.7125G>T'),
     ('NG_012337.1(NM_003002.2):c.1del',
      'NG_012337.1:g.5062del'),
     ('NG_012337.1(NM_003002.2):c.-1del',
      'NG_012337.1:g.5061del'),
     ('NG_012337.1(NM_003002.2):c.52+1del',
      'NG_012337.1:g.5114del'),
     ('NG_012337.1(NM_003002.2):c.*824del',
      'NG_012337.1:g.13948del'),
     ('NG_012337.1(NM_003002.2):c.*824+10del',
      'NG_012337.1:g.13958del'),

     # To be curated.
     # --------------
     # - Fix the ambiguity for the length/location.
     # ('NG_012337.1:g.100_200>400',
     #  'NG_012337.1:g.100_200delins100'),
     ])
def test_mutalyzer3(hgvs_description, normalized_description, monkeypatch):
    monkeypatch.setattr('retriever.retriever.fetch_annotations',
                        fetch_annotation)
    monkeypatch.setattr('retriever.retriever.fetch_sequence',
                        fetch_sequence)
    assert mutalyzer3(hgvs_description) == normalized_description
