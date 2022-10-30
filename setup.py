from setuptools import setup

setup(
    name="gmxpla",
    version="0.0.1",
    install_requires=[
        "pyyaml", "matplotlib",
    ],
    entry_points={
        'console_scripts': [
            'gmxpla=gmxpla.gmxpla:main',
        ],
    },
    author="Michio Katouda",
    author_email="katouda@rist.or.jp",
    description="gromax protein-ligand MD trajectory analysis tools",
    url="https://github.com/mkatouda/gmxpla",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.7',
)
