# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

include_directories(${SOSIM_INCLUDE_DIR})

add_subdirectory(Framework)

set(sosim-shared-libs "")
add_subdirectory(Shared)

add_subdirectory(DFSPH)

add_subdirectory(MSMSPH)

#list(LENGTH sosim-shared-libs length)
#message(${length})