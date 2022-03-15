import DS

VERSION = DS.__version__
DISTNAME = "Data science for environmental sciences"
DESCRIPTION = ("Module for processing arrays of simulation results, analysing them, generating graphical outputs and solving reaction networks",)
AUTHORS = "Swamini Khurana"
EMAIL = "swamini.khurana@gmail.com"
LICENSE = "Copyright(c) 2005-2020, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. All rights reserved."
PROJECT_URLS = {
    "Github": "https://github.com/Skhurana136/dsenvsci",
    "Documentation": "TBDone",
}


from setuptools import setup, find_packages

setup(
    name="dsenvsci",
    version="0.0",
    description='sk_package_for_data_science',
    authors='Swamini Khurana',
    author_email='swamini.khurana@gmail.com',
    packages=find_packages(exclude=["tests*", "docs*"]),
    license="Copyright(c) 2005-2020, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. All rights reserved.",
)
