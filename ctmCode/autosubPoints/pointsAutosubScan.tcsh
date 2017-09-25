#!/bin/tcsh

foreach i(`seq 80 120`)
  cp autosubPoints.txt $i.txt
  sed -i 's/X/'$i'/g' $i.txt
end
