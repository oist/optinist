from setuptools import setup

def _requires_from_file(filename):
    return open(filename).read().splitlines()

setup(
    install_requires=_requires_from_file('requirements.txt'),
    tests_require=["pytest"],
    # package_dir={"": "src"},
    entry_points = {
        'console_scripts': [
            'run_optinist=main:main',
        ]
    },
    extras_require = {
        "docs": [
        'sphinx>=4.3.0,<4.6.0',
        'sphinxcontrib-apidoc',
        'sphinx_rtd_theme',
        'sphinx-prompt',
        'sphinx-autodoc-typehints', # Automatically adds types to docs
        'sphinx-copybutton==0.5.0',
        'sphinx-autobuild==2021.3.14' # Live rebuilding and reloading of docs for developing locally.
        ],
    },
)
