from setuptools import setup, find_packages

setup(
    name="asodesigner",
    version="0.1.0",
    author="TAU_iGEM_2025",
    author_email="michaelkovaliov97@gmail.com",
    description="A short description of asodesigner",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/RedPenguin100/ASOdesign",
    packages=find_packages(),  # <- automatically finds asodesigner/
    install_requires=[
        "numpy",
        "viennarna",
        "pandas",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
