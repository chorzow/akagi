missed = $(shell flakeheaven missed)

before-lint =
# if not all linter plugins are installed (flakehell missed returns a non-empty string),
# add linter-plugins to lint prerequisites
ifneq ($(missed),)
    before-lint = linter-plugins
endif

lint: $(before-lint)
	flakeheaven lint components

linter-plugins:
	pip install $(missed) -r requirements-dev.txt

requirements-dev:
	pip install -r requirements-dev.txt