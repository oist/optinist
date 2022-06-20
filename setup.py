from setuptools import setup

def _requires_from_file(filename):
    return open(filename).read().splitlines()

setup(
    install_requires=_requires_from_file('requirements.txt'),
    tests_require=["pytest"],
    entry_points = {
        'console_scripts': [
            'run_optinist=main:main',
        ]
    },
)
