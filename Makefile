.PHONY: test_run
test_run:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build test_studio
	docker-compose -f docker-compose.yml build test_studio_frontend
	docker-compose -f docker-compose.yml run test_studio
	docker-compose -f docker-compose.yml run test_studio_frontend

.PHONY: test_python
test_python:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build test_studio
	docker-compose -f docker-compose.yml run test_studio

.PHONY: test_frontend
test_frontend:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build test_studio_frontend
	docker-compose -f docker-compose.yml run test_studio_frontend

.PHONY: build_frontend
build_frontend:
	docker-compose -f docker-compose.yml down --rmi all --volumes --remove-orphans
	docker-compose -f docker-compose.yml rm -f
	docker-compose -f docker-compose.yml build build_studio_frontend
	docker-compose -f docker-compose.yml run build_studio_frontend

.PHONY: docs
docs:
	rm -rf docs/_build/
	# sphinx-apidoc -f -o ./docs/_build/modules ./studio
	sphinx-autobuild -b html docs docs/_build --port 8001

.PHONY: dockerhub
dockerhub:
	docker build --rm -t oistncu/studio:latest . --platform=linux/amd64
	docker push oistncu/studio:latest

.PHONY: local_build
local_build:
	cp -r frontend/build studio/frontend/build
	pip install .

.PHONY: upload_testpypi
upload_testpypi:
	python setup.py sdist
	twine upload --repository testpypi dist/*

.PHONY: test_pypi
test_pypi:
	python3 -m pip install --index-url https://test.pypi.org/simple/ studio

.PHONY: push_pypi
push_pypi:
	twine upload --repository pypi dist/*