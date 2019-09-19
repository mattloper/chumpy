all:

upload:
	python setup.py sdist
	twine upload dist/*


#sdist:
#	python setup.py sdist && rsync -avz dist/chumpy-0.5.tar.gz files:~/chumpy/latest.tgz && python ./api_compatibility.py && rsync -avz ./api_compatibility.html files:~/chumpy/

clean: 

tidy: 

test:  clean qtest
qtest:   all
	# For some reason the import changes for Python 3 caused the Python 2 test
	# loader to give up without loading any tests. So we discover them ourselves.
	# python -m unittest
	find chumpy -name 'test_*.py' | sed -e 's/\.py$$//' -e 's/\//./' | xargs python -m unittest

coverage: clean qcov
qcov: all
	env LD_PRELOAD=$(PRELOADED) coverage run --source=. -m unittest discover -s .
	coverage html
	coverage report -m

