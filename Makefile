.PHONY: test_run
test_run:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build
	docker-compose -f docker-compose.test.yml run optinist
	docker-compose -f docker-compose.test.yml run optinist_frontend

.PHONY: test_python
test_python:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build
	docker-compose -f docker-compose.test.yml run optinist

.PHONY: test_frontend
test_frontend:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build
	docker-compose -f docker-compose.test.yml run optinist_frontend

.PHONY: docs
docs:
	rm -rf docs/_build/
	# sphinx-apidoc -f -o ./docs/_build/modules ./optinist
	sphinx-autobuild -b html docs docs/_build --port 8001
