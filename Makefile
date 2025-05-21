all:

upload:
	uv publish

test:
	uv run python -m unittest

coverage: clean qcov
qcov: all
	uv run coverage run -m unittest
	uv run coverage html
	uv run coverage report -m
