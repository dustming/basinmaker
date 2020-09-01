import pytest

def pytest_addoption(parser):
    parser.addoption("--HYDAT_Path", action="store")
    parser.addoption("--name", action="store", default="default name")



@pytest.fixture(scope='session')
def HYDAT_Path(request):
    name_value = request.config.option.HYDAT_Path
    if name_value is None:
        pytest.skip()
    return name_value    