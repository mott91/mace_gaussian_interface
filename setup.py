from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = []
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Skip editable or git dependencies
        if line.startswith("-e") or re.match(r"git\+https?://", line):
            continue
        requirements.append(line)

setup(
    name="mace-gaussian-interface",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@domain.com",
    description="Interface between MACE ML potentials and Gaussian for molecular calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/your-repo-name",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "gm-main=gm_main:main",
        ],
    },
)
