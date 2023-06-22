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
	docker-compose -f docker-compose.yml build test_optinist
	docker-compose -f docker-compose.yml run test_optinist

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
	pip install -e '.[doc]'
	# sphinx-apidoc -f -o ./docs/_build/modules ./optinist
	sphinx-autobuild -b html docs docs/_build --port 8001

.PHONY: dockerhub
dockerhub:
	docker build --rm -t oistncu/optinist:latest . --platform=linux/amd64
	docker push oistncu/optinist:latest

.PHONY: local_build
local_build:
	cd frontend
	yarn install && yarn build
	cd ../
	pip install .

.PHONY: upload_testpypi
upload_testpypi:
	python -m build
	twine upload --repository testpypi dist/*

.PHONY: test_pypi
test_pypi:
	python3 -m pip install --index-url https://test.pypi.org/simple/ optinist

.PHONY: push_pypi
push_pypi:
	twine upload --repository pypi dist/*


.PHONY: format
format:
	black optinist *.py
	isort optinist *.py
	flake8 optinist *.py
