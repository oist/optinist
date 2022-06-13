.PHONY: test_run
test_run:
	python -m pytest optinist/
	cd frontend && rm -rf node_modules && yarn install && yarn test:ci

.PHONY: test_python
test_python:
	python -m pytest optinist/


.PHONY: docs
docs:
	rm -rf docs/_build/
	# sphinx-apidoc -f -o ./docs/_build/modules ./optinist
	sphinx-autobuild -b html docs docs/_build --port 8001
