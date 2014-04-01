all:


sdist:
	python setup.py sdist && rsync -avz dist/chumpy-0.5.tar.gz files:~/chumpy/latest.tgz && python ./api_compatibility.py && rsync -avz ./api_compatibility.html files:~/chumpy/

clean: 

tidy: 

test:  clean qtest
qtest:   all
	python -m unittest discover

coverage: clean qcov
qcov: all
	env LD_PRELOAD=$(PRELOADED) coverage run --source=. -m unittest discover -s .
	coverage html
	coverage report -m

