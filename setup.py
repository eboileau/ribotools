#! /usr/bin/env python3

from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop


class SetupInstall(install):

    def run(self):
        install.run(self)


class SetupDevelop(develop):

    def run(self):
        develop.run(self)


setup(
    cmdclass={
        'install': SetupInstall,
        'develop': SetupDevelop
    }
)
