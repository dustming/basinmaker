import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="basinmaker", # Replace with your own username
    version="0.0.1",
    author="Ming Han",
    author_email="m43han@uwaterloo.ca",
    description="A GIS toolbox to delineate watershed with lakes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dustming/basinmaker",
    packages = setuptools.find_packages(
        where = 'basinmaker',
    ),
    package_dir = {"":"basinmaker"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)