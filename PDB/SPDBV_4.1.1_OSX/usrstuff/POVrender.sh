#!/bin/sh
if [ $# -gt 5 ]; then
if [ $6 == "-s" ]; then

# --------- STEREO ----------------------------------

# retrieve the line number of the rotate camera commands

	leftnum=`grep -n '// rotate' $5 | head -1 | awk -F: '{print $1}'`
	rightnum=`grep -n '// rotate' $5 | tail -1 | awk -F: '{print $1}'`

# uncomment the camera rotation for left and right eyes

	sed -e "s/right <100\/100,/right <$2\/$4,/" -e "$leftnum s/\/\/ rotate /rotate /1" $5 > /tmp/left$$.pov
	sed -e "s/right <100\/100,/right <$2\/$4,/" -e "$rightnum s/\/\/ rotate /rotate /1" $5 > /tmp/right$$.pov

#render each view

	/usr/local/bin/povray +W$2 +H$4 +I/tmp/left$$.pov +O/tmp/left$$.png
	/usr/local/bin/povray +W$2 +H$4 +I/tmp/right$$.pov +O/tmp/right$$.png

# open files

	mv /tmp/left$$.png $5left.png
	mv /tmp/right$$.png $5right.png
	open $5left.png $5right.png

exit;
fi
fi

# --------- NORMAL ----------------------------------

sed -e "s/right <100\/100,/right <$2\/$4,/" $5 > /tmp/$$.pov
/usr/local/bin/povray +W$2 +H$4 +I/tmp/$$.pov +O/tmp/$$.png
mv /tmp/$$.png $5.png
open $5.png
