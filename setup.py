import setuptools

long_message = 'Fast Response Analysis'
version = "0.0.1"

setuptools.setup(
    name="fast_response", 
    version=version,
    author="Pizzuto, Alex",
    author_email="",
    description="Code for performing rapid neutrino followup",
    long_description=long_message,
    #long_description_content_type="text/markdown",
    url="https://github.com/icecube/FastResponseAnalysis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
    ],
    python_requires='>=3.1',
)
