#TODO : change minimum version?
cmake_minimum_required(VERSION 3.6.0)
project(oxb1tz)



include_directories(
#sources
	${PROJECT_SOURCE_DIR}/include
	${PROJECT_SOURCE_DIR}/src
#unit test
	${PROJECT_SOURCE_DIR}/unit_tests
)


file(GLOB SRCS
#sources
	"${PROJECT_SOURCE_DIR}/include/*"
	"${PROJECT_SOURCE_DIR}/src/*"
	"${PROJECT_SOURCE_DIR}/src/oxb1tz/*"
	"${PROJECT_SOURCE_DIR}/src/oxb1tz/core2d/*"

#unit test
	"${PROJECT_SOURCE_DIR}/unit_tests/*"
)

#add libs
add_library(oxb1tz
 ${SRCS}
)

#compiler optimization flags
target_compile_options(oxb1tz PRIVATE
#release
 $<$<COMPILE_LANGUAGE:C>: -DOXB_BUILD -Wall -ffast-math -Wextra -Werror -O3 -fvectorize -std=c11>
#just default, no optimization
# $<$<COMPILE_LANGUAGE:C>: -DOXB_BUILD -O0>
)
