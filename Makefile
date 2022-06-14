.PHONY: test_run
test_run:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build
	docker-compose -f docker-compose.test.yml up
	# cd frontend && rm -rf node_modules && yarn install && yarn test:ci

.PHONY: test_python
test_python:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build
	docker-compose -f docker-compose.test.yml up

.PHONY: docs
docs:
	rm -rf docs/_build/
	# sphinx-apidoc -f -o ./docs/_build/modules ./optinist
	sphinx-autobuild -b html docs docs/_build --port 8001
