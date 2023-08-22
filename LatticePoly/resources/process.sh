##
##  process.sh
##  LatticePoly
##
##  Created by mtortora on 27/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

#!/bin/bash -l

if [ "$#" -ge "3" ]; then
	find $2 -type d -links 2 -exec python3 $1 {} ${@: 3} \;
else
	echo -e "\033[1;31mUsage is $0 scriptName rootDataDir initFrame [opt]\033[0m"
fi
