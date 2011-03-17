.SUFFIXES:
.PHONY: FORCE
.DELETE_ON_ERROR:

default:
	@echo "ain't no $@ target" 1>&2; exit 1


setup:
	ln -s /usr/lib/pymodules/python2.6/Bio .


.PHONY: clean cleaner cleanest
clean:
	find . \( -name '*~' -o -name '*.bak' \) -print0 | xargs -0r rm -f

cleaner: clean
	find . -name '*.pyc' -print0 | xargs -0r rm -f

cleanest: cleaner
