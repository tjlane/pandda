import pytest

def test_DatasetProcessor():

    def f1(d):
        d.add(1)
        return d
    def f2(d):
        d.add(2)
        return d
    def f3(d):
        d.add(3)
        return d

    from giant.mulch.utils import DatasetProcessor

    processor = DatasetProcessor(
        functions = [f1,f2,f3],
        )

    s = set()
    s_out = processor(s)

    assert (s is s_out)
    assert sorted(list(s)) == [1,2,3]
