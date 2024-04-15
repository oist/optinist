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

.PHONY: format
format:
	black studio *.py
	isort studio *.py
	flake8 studio *.py

.PHONY: docs
docs:
	rm -rf docs/_build/
	poetry install --with doc --no-root
	# sphinx-apidoc -f -o ./docs/_build/modules ./studio
	sphinx-autobuild -b html docs docs/_build --port 8001

.PHONY: local_build
local_build:
	rm -rf dist
	cd frontend && yarn install --ignore-scripts && yarn build
	poetry build

.PHONY: push_testpypi
push_testpypi:
	poetry publish -r testpypi

.PHONY: install_testpypi
install_testpypi:
	pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ optinist==${ver}
	pip show optinist

.PHONY: build_test_docker
build_test_docker:
	docker build --no-cache -t optinist-release-test:${ver} -f studio/config/docker/Dockerfile .

.PHONY: run_test_docker
run_test_docker:
	docker run -it \
	-v ${volume}:/app/studio_data \
	--env OPTINIST_DIR="/app/studio_data" \
	--name optinist-release-test -d -p 8000:8000 optinist-release-test:${ver} \
	poetry run python main.py --host 0.0.0.0 --port 8000

.PHONY: push_pypi
push_pypi:
	poetry publish

.PHONY: push_dockerhub
push_dockerhub:
	docker build --rm -t oistncu/optinist:latest -f studio/config/docker/Dockerfile . --platform=linux/amd64
	docker push oistncu/optinist:latest
