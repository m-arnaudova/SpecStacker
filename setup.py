# setup.py

from setuptools import setup, find_packages

# Read the contents of your README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='SpecStacker',  # Replace with your own package name
    version='1.0.0',  # Initial release version
    author='Marina Arnaudova',  
    author_email='m.i.arnaudova@gmail.com',  
    description='This is a new rest-frame spectral stacking code.',  
    long_description=long_description,  # Long description read from the README file
    url='https://github.com/m-arnaudova/SpecStacker', 
    packages=find_packages(),  # Automatically find packages in the current directory
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Specify the Python versions you support
    install_requires=[
        # List your package dependencies here, e.g.,
        'astropy','extinction','lmfit','matplotlib','numpy','scipy','sfdmap2','spectres'
    ],
)
