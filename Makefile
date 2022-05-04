.PHONY: test_run
test_run:
	python -m pytest
	cd frontend && rm -rf node_modules && yarn install && yarn test:ci