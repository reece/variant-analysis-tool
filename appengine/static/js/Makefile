.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

default:
	@echo "no $@ target" 1>&2; exit 1

update: jquery jquery-ui

jquery: %: /usr/share/javascript/%
	/bin/cp -abv $< .

jquery-ui: %: /usr/share/javascript/%
	/bin/cp -abv $< .
