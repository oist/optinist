.PHONY: test_run
test_run:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build test_studio
	docker-compose -f docker-compose.test.yml build test_studio_frontend
	docker-compose -f docker-compose.test.yml run test_studio
	docker-compose -f docker-compose.test.yml run test_studio_frontend

.PHONY: test_python
test_python:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build test_studio
	docker-compose -f docker-compose.test.yml run test_studio

.PHONY: test_frontend
test_frontend:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build test_studio_frontend
	docker-compose -f docker-compose.test.yml run test_studio_frontend

.PHONY: build_frontend
build_frontend:
	docker-compose -f docker-compose.test.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.test.yml rm -f
	docker-compose -f docker-compose.test.yml build build_studio_frontend
	docker-compose -f docker-compose.test.yml run build_studio_frontend

.PHONY: docs
docs:
	rm -rf docs/_build/
	poetry install --with doc --no-root
	# sphinx-apidoc -f -o ./docs/_build/modules ./studio
	sphinx-autobuild -b html docs docs/_build --port 8001

.PHONY: dockerhub
dockerhub:
	docker build --rm -t oistncu/optinist:latest . --platform=linux/amd64
	docker push oistncu/optinist:latest

.PHONY: local_build
local_build:
	cd frontend && yarn install --ignore-scripts && yarn build
	poetry build

.PHONY: upload_testpypi
upload_testpypi:
	poetry publish -r testpypi

.PHONY: install_testpypi
install_testpypi:
	pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ optinist
	pip list | grep optinist

.PHONY: push_pypi
push_pypi:
	poetry publish

.PHONY: format
format:
	black studio *.py
	isort studio *.py
	flake8 studio *.py
