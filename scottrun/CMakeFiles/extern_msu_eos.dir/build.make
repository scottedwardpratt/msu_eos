# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.26.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.26.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/scottpratt/git/msu_eos/scottrun

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/scottpratt/git/msu_eos/scottrun

# Utility rule file for extern_msu_eos.

# Include any custom commands dependencies for this target.
include CMakeFiles/extern_msu_eos.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/extern_msu_eos.dir/progress.make

CMakeFiles/extern_msu_eos:
	cd /Users/scottpratt/git/msu_eos/software && make

extern_msu_eos: CMakeFiles/extern_msu_eos
extern_msu_eos: CMakeFiles/extern_msu_eos.dir/build.make
.PHONY : extern_msu_eos

# Rule to build all files generated by this target.
CMakeFiles/extern_msu_eos.dir/build: extern_msu_eos
.PHONY : CMakeFiles/extern_msu_eos.dir/build

CMakeFiles/extern_msu_eos.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/extern_msu_eos.dir/cmake_clean.cmake
.PHONY : CMakeFiles/extern_msu_eos.dir/clean

CMakeFiles/extern_msu_eos.dir/depend:
	cd /Users/scottpratt/git/msu_eos/scottrun && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/scottpratt/git/msu_eos/scottrun /Users/scottpratt/git/msu_eos/scottrun /Users/scottpratt/git/msu_eos/scottrun /Users/scottpratt/git/msu_eos/scottrun /Users/scottpratt/git/msu_eos/scottrun/CMakeFiles/extern_msu_eos.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/extern_msu_eos.dir/depend

