# set up environment variables
# the ${VAR:+:} part adds a double colon only if VAR is not empty
export PATH="/home/lreichen/work/acts/acts/install/bin${PATH:+:}${PATH}"
export LD_LIBRARY_PATH="/home/lreichen/work/acts/acts/install/lib64${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="/home/lreichen/work/acts/acts/install/lib64${DYLD_LIBRARY_PATH:+:}${DYLD_LIBRARY_PATH}"
export CMAKE_PREFIX_PATH="/home/lreichen/work/acts/acts/install${CMAKE_PREFIX_PATH:+:}${CMAKE_PREFIX_PATH}"
