.PHONY: test_run
test_run:
	python -m pytest
	cd frontend && rm -rf node_modules && yarn install && yarn test:ci

.PHONY : docs
docs:
	rm -rf docs/build/
	# sphinx-apidoc -f -o ./docs/_build/modules ./optinist
	sphinx-autobuild -b html docs docs/_build --port 8001