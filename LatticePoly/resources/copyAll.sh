##
##  copyAll.sh
##  LatticePoly
##
##  Created by mtortora on 02/04/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

#!/bin/bash -l

if [ "$#" -eq "2" ]; then
	find $1 -type f ! \( -iname "*.vtp" -o -iname "*.tar.gz" \) -exec bash -c 'mkdir -p "$0"/$(dirname {}) && cp {} "$0"/{}' $2 \;
else
	echo -e "\033[1;31mUsage is $0 rootDataDir copyDataDir\033[0m"
fi
