.PHONY: test_run
test_run:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build test_optinist
	docker-compose -f docker-compose.yml build test_optinist_frontend
	docker-compose -f docker-compose.yml run test_optinist
	docker-compose -f docker-compose.yml run test_optinist_frontend

.PHONY: test_python
test_python:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build optinist
	docker-compose -f docker-compose.yml run optinist

.PHONY: test_frontend
test_frontend:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build test_optinist_frontend
	docker-compose -f docker-compose.yml run test_optinist_frontend

.PHONY: build_frontend
build_frontend:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build build_optinist_frontend
	docker-compose -f docker-compose.yml run build_optinist_frontend

.PHONY: docs
docs:
	rm -rf docs/_build/
	# sphinx-apidoc -f -o ./docs/_build/modules ./optinist
	sphinx-autobuild -b html docs docs/_build --port 8001
