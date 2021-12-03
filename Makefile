SHELL := /bin/bash

# This is going to be used in deployment
APP := "sxdm_py38"

TSTAMP:=$(shell date +"%Y%m%d-%H%M")
TF:="/sls/MX/applications/$(APP)"
SOURCE:="./src/"

UIFILES := $(wildcard *.ui)
PYFILES := $(UIFILES:.ui=.py)

TBL := $(shell echo $(BEAMLINE_XNAME) | tr A-Z a-z)

ECHO = /bin/echo

.PHONY: all-beamlines

check:  ui git git-dirty git-tagged

all-beamlines: git git-tagged
		$(eval B := $@)
		$(eval T:=$(TF)/$(B)/$(TARGET))
		$(eval STABLE:=$(TF)/$(B)/stable)
		@echo "Deploying app=$(APP) as $(TARGET) for $(B) on folder: $(T)"
		@if test -L $(STABLE) ; then /bin/rm $(STABLE) ; fi
		rsync -av --exclude 'tester*' \
						--exclude '*~' --exclude '*pyc' \
						--exclude Makefile \
						--exclude .git* \
						--exclude .idea \
						--exclude README.md \
						--delete --delete-excluded $(SOURCE) $(T)
		@cd $(TF)/$(B) && /bin/ln -s $(TARGET) stable
		@echo Latest stable is now: $(T)
		@echo "...finished."

git-bogus:
		$(eval TARGET := devel)

git:
		$(eval COMMIT_TAG := $(shell git describe --exact-match HEAD))
		$(eval COMMIT_ID := $(shell git rev-parse --short HEAD))
		$(eval COMMIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD))
		$(eval TARGET:=$(COMMIT_TAG)-$(COMMIT_ID)-$(COMMIT_BRANCH))

git-dirty: ## (internal) used to check if the git repository is dirty
		@echo Checking git status
		@status=$$(git status --porcelain); \
		if test "x$${status}" != x; then \
				echo "Working directory is dirty -- commit and tag before deploying" >&2; \
				exit 2 ; \
		fi

git-tagged: ## (internal) used to check if the last commit is also annotated (MANDATORY)
		@echo Checking git tag
		git describe --exact-match HEAD 2>/dev/null || ( echo "Code is not tagged -- tag before deploying"; exit 123 )

amend-git: ## amend previous commit keeping message and add everything modified so far
		git commit -a --amend --no-edit

tag: ##git-dirty ## use to make a standard TAG, repository must be clean
		git tag -f -a -m "Version $(TSTAMP)" v$(TSTAMP)
		#git push --tags

mtag: git-dirty ## use to make a standard TAG, repository must be clean
		git tag -f -a v$(TSTAMP)

ui: $(PYFILES)

%.py: %.ui
		pyuic5 $< --output $@

