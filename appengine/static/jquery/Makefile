.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

V:=1.8.2

default:
	@echo "no $@ target" 1>&2; exit 1

.PHONY: update
update: update-$V.log

update-$V.log: /tmp/jquery-ui-$V.custom.zip
	set -x; unzip -x $< css/smoothness/\* js/\* 1>$@ 2>&1

/tmp/jquery-ui-$V.custom.zip:
	wget -P /tmp/ http://jqueryui.com/download/jquery-ui-$V.custom.zip

