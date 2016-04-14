from distutils.core import setup

setup(
    name='occult',
    version='0.1dev',
    packages=['occult','occult.model','occult.observing'],
    license='MIT License',
    long_description=open('README.md').read(),
)