#!/bin/bash
# echo "python -m pip install --upgrade pip"
# python -m pip install --upgrade pip
# echo "python -m pip install --upgrade build twine"
# python -m pip install --upgrade build twine
echo
echo "Cleaning up dist directory."
echo
rm -rf dist/*
echo "python -m build"
python -m build
echo "python -m twine upload --repository testpypi dist/*"
python -m twine upload --repository testpypi dist/*
echo
echo "Test install:"
echo "python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps scCirRL"
echo
echo "All good? Upload to real PyPI:"
echo python3 -m twine upload dist/*
echo
