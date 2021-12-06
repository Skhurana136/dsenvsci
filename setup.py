import DS

VERSION = DS.__version__
DISTNAME = "Data science for environmental sciences"
DESCRIPTION = ("Module for processing arrays of simulation results, analysing them, generating graphical outputs and solving reaction networks",)
AUTHORS = "Swamini Khurana"
EMAIL = "swamini.khurana@gmail.com"
LICENSE = "Copyright(c) 2005-2020, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. All rights reserved."
PROJECT_URLS = {
    "GitLab": "https://git.ufz.de/ml-cafe/ml-cafe_project_soilmoisture/-/tree/master/SM_module",
    "Documentation": "TBDone",
}


from setuptools import setup, find_packages

setup(
    name=DISTNAME,
    version=VERSION,
    description=DESCRIPTION,
    authors=AUTHORS,
    author_email=EMAIL,
    packages=find_packages(exclude=["tests*", "docs*"]),
    license=LICENSE,
)
