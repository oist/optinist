import os

from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, "optinist", "version.py")).read())


def _requires_from_file(filepath):
    def take_package_name(name):
        # if "snakemake" in name:
        #     return (
        #         "snakemake @ https://github.com/ShogoAkiyama/snakemake/"
        #         "archive/refs/tags/v7.7.2-post1.zip"
        #     )
        # else:
        return name.strip()

    with open(filepath) as fp:
        return [take_package_name(pkg_name) for pkg_name in fp.readlines()]


setup(
    name="optinist",
    version=VERSION,  # noqa: F821
    description="Calcium Imaging Pipeline Tool",
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
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS :: MacOS X",
    ],
    install_requires=_requires_from_file("requirements.txt"),
    packages=find_packages(exclude=["optinist/tests*"]),
    tests_require=["pytest"],
    zip_safe=False,
    include_package_data=True,
    package_data={
        "": [
            "frontend/build/*",
            "frontend/build/static/*",
            "frontend/build/static/css/*",
            "frontend/build/static/js/*",
            "frontend/build/static/media/*",
            "app_config/*.yaml",
            "config/*.yaml",
            "conda/*.yaml",
            "Snakefile",
        ]
    },
    py_modules=["Snakefile"],
    entry_points={
        "console_scripts": [
            "run_optinist=optinist.__main__:main",
        ]
    },
)
