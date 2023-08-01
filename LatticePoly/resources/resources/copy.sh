##
##  copy.sh
##  LatticePoly
##
##  Created by mtortora on 02/04/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

#!/bin/bash -l

if [ "$#" -eq "3" ]; then
	find $2 -iname $1 -exec bash -c 'mkdir -p "$0"/$(dirname {}) && cp {} "$0"/{}' $3 \;
else
	echo -e "\033[1;31mUsage is $0 fileName rootDataDir copyDataDir\033[0m"
fi
