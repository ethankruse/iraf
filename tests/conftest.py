import pytest


@pytest.fixture(scope='session')
def combine_dir(tmpdir_factory):
    outdir = tmpdir_factory.mktemp('combine')
    return outdir
