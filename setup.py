import setuptools

long_message = 'Fast Response Analysis'
version = "1.0.0"

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
    install_requires=[
        'astropy==2.0.16',
        'healpy==1.13.0',
        'matplotlib==2.2.5',
        'numpy==1.16.6',
        'pandas==0.24.2',
        'pyfiglet==0.8.post1',
        'python-dateutil==2.8.1',
        'pyzmq==19.0.1',
        'scipy==1.2.3',
        'seaborn==0.9.1',
        'zmq==0.0.0',
        'py27hash==1.0.2',
        'pygcn==1.0.2'
    ]
)
