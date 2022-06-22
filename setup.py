from setuptools import setup, find_packages, Extension

def _requires_from_file(filepath):

    def take_package_name(name):
        if name.startswith("git+"):
            _package = name.strip("=")[-1]
            return f"{_package}"
        else:
            return name.strip()

    with open(filepath) as fp:
        return [take_package_name(pkg_name) for pkg_name in fp.readlines()]

setup(
    name="optinist",
    version="0.1.0-6",
    description="An offline deep reinforcement learning library",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/oist/optinist",
    author="OIST",
    license="GPL3.0",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: POSIX :: Linux",
        'Operating System :: Microsoft :: Windows',
        "Operating System :: MacOS :: MacOS X"
    ],
    install_requires=_requires_from_file('requirements.txt'),
    dependency_links=["git+https://github.com/ShogoAkiyama/snakemake@main#egg=snakemake"],
    packages=find_packages(exclude=["optinist/tests*"]),
    tests_require=["pytest"],
    zip_safe=False,
    include_package_data=True,
    package_data={"": [
        "frontend/build/*",
        "frontend/build/static/*",
        "frontend/build/static/css/*",
        "frontend/build/static/js/*",
        "frontend/build/static/media/*",
        "config/*.yaml",
        "conda/*.yaml",
    ]},
    py_modules=["Snakefile"],
    entry_points={
        "console_scripts": [
            "run_optinist=optinist.__main__:main",
        ]
    },
)
