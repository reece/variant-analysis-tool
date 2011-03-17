.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

LOGO_DIR:=${HOME}/Pictures/logos
FMT:=png

.PHONY: all
all: qb3-logo.${FMT} celera-logo.${FMT} berkeley-logo.${FMT} ccb-text-logo.${FMT}
	identify -format '<img src="/static/gfx/%f" width=%w height=%h target="_blank" />\n' $^

qb3-logo.${FMT}: ${LOGO_DIR}/qb3.png
celera-logo.${FMT}: ${LOGO_DIR}/celera.png
berkeley-logo.${FMT}: ${LOGO_DIR}/berkeley.png
ccb-text-logo.${FMT}: ${LOGO_DIR}/ccb-text.png


# convert to 30px tall logo (logo will be as wide as needed, max 200)
%.${FMT}:
	convert -geometry 200x30 $< $@

