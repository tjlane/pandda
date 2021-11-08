import os, copy
import pytest
import pathlib as pl

from giant.tests.jiffies.resources import (
    MultiStateResources,
    )

@pytest.fixture(scope="module")
def resources():

    this = pl.Path(__file__)

    resources_path = (
        this.parent / 'resources'
        )

    assert resources_path.exists()

    reso = MultiStateResources(resources_path)

    return reso
