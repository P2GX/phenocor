# Internal

This page contains internal instructions for maintaining and releasing the project.



## Documentation

These pages are generated with Sphinx.

To set things up, perform the following steps (substitute name of venv if needed).

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -e ".[docs]"
```

Then build the docs:
```bash
cd docs
make html
```

Open:
docs/build/html/index.html
 

# PyPI release

Install the packages build and twine. Then

```bash
python -m build
python -m twine upload dist/*
```