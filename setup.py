import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="basinmaker", 
    version="2.0.6",
    author="basinmaker development team",
    author_email="m43han@uwaterloo.ca",
    description="An automated GIS toolbox for watershed delineation with lakes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dustming/basinmaker",
    packages = setuptools.find_packages(
        include=[
            "basinmaker",
            "basinmaker.*",
        ],
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)