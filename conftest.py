try:
    import giant.logs as lg
    # uncomment to enable logging in all modules during testing
    #logger = lg.setup_logging_basic(__name__)
except Exception as e:
    pass

import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--fast", action="store_true", default=False, help="run only fast tests"
    )

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")

def pytest_collection_modifyitems(config, items):
    if config.getoption("--fast"):
        skip_slow = pytest.mark.skip(reason="will not run with --fast option")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)
