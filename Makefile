
DIR='/exchange/mx/daplus/sxdm'

PREFIX=$1
SXDM=$(PREFIX)

sls:
	test -d $(DIR) || mkdir $(DIR);
	cd src && rsync -av --exclude '*.ipynb*' --exclude '*~' --exclude '*pyc' --exclude 'test.py' --delete --delete-excluded . $(DIR)

diff:
	cd src && for f in `find . -name '*.py'`; do echo -e "### diff $$f $(DIR)/$$f"; diff $$f $(DIR)/$$f 2> /dev/null ; done

install:
	test -d $(SXDM) || mkdir -p $(SXDM);
	test -d $(SXDM)/lib || mkdir $(SXDM)/lib;
	cd src && rsync -av --exclude '*.ipynb*' --exclude '*~' --exclude '*pyc' --exclude 'test.py' --delete --delete-excluded . $(SXDM)/lib/
	cd .. 
	rsync -av guis $(SXDM)/.
	rsync -av bin $(SXDM)/.
	rsync -av xds $(SXDM)/.
