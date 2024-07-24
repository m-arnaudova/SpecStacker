# setup.py

from setuptools import setup, find_packages

# Read the contents of your README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='SpecStacker',  # Replace with your own package name
    version='0.1.0',  # Initial release version
    author='Marina Arnaudova',  # Replace with your name
    author_email='m.i.arnaudova@example.com',  # Replace with your email
    description='This is a new rest-frame spectral stacking code.',  # Short description of your package
    long_description=long_description,  # Long description read from the README file
    long_description_content_type='text/markdown',  # Optional (see note above)
    url='https://github.com/m-arnaudova/SpecStacker',  # Replace with the URL of your package
    packages=find_packages(),  # Automatically find packages in the current directory
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Specify the Python versions you support
    install_requires=[
        # List your package dependencies here, e.g.,
        'astropy','extinction','lmfit','matplotlib','numpy','scipy','sfdmap','spectres'
    ],
)
