##
##  unpack.sh
##  LatticePoly
##
##  Created by mtortora on 27/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

#!/bin/bash -l

if [ "$#" -eq "1" ]; then
	find $1 -iname "*.tar.xz" -execdir tar -xzvmf {} \;
else
	echo -e "\033[1;31mUsage is $0 rootDataDir\033[0m"
fi
