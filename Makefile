.PHONY: test_run
test_run:
	python -m pytest
	cd frontend && yarn test