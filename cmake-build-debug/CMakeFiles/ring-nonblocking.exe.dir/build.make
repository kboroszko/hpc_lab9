# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/kajetan/IntellijCLion/clion-2019.2.5/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/kajetan/IntellijCLion/clion-2019.2.5/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kajetan/studia/hpc/hpc_lab9

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kajetan/studia/hpc/hpc_lab9/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/ring-nonblocking.exe.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ring-nonblocking.exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ring-nonblocking.exe.dir/flags.make

CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.o: CMakeFiles/ring-nonblocking.exe.dir/flags.make
CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.o: ../ring-nonblocking.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kajetan/studia/hpc/hpc_lab9/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.o"
	CC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.o -c /home/kajetan/studia/hpc/hpc_lab9/ring-nonblocking.cpp

CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.i"
	CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kajetan/studia/hpc/hpc_lab9/ring-nonblocking.cpp > CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.i

CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.s"
	CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kajetan/studia/hpc/hpc_lab9/ring-nonblocking.cpp -o CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.s

# Object files for target ring-nonblocking.exe
ring__nonblocking_exe_OBJECTS = \
"CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.o"

# External object files for target ring-nonblocking.exe
ring__nonblocking_exe_EXTERNAL_OBJECTS =

ring-nonblocking.exe: CMakeFiles/ring-nonblocking.exe.dir/ring-nonblocking.cpp.o
ring-nonblocking.exe: CMakeFiles/ring-nonblocking.exe.dir/build.make
ring-nonblocking.exe: CMakeFiles/ring-nonblocking.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/kajetan/studia/hpc/hpc_lab9/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ring-nonblocking.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ring-nonblocking.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ring-nonblocking.exe.dir/build: ring-nonblocking.exe

.PHONY : CMakeFiles/ring-nonblocking.exe.dir/build

CMakeFiles/ring-nonblocking.exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ring-nonblocking.exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ring-nonblocking.exe.dir/clean

CMakeFiles/ring-nonblocking.exe.dir/depend:
	cd /home/kajetan/studia/hpc/hpc_lab9/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kajetan/studia/hpc/hpc_lab9 /home/kajetan/studia/hpc/hpc_lab9 /home/kajetan/studia/hpc/hpc_lab9/cmake-build-debug /home/kajetan/studia/hpc/hpc_lab9/cmake-build-debug /home/kajetan/studia/hpc/hpc_lab9/cmake-build-debug/CMakeFiles/ring-nonblocking.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ring-nonblocking.exe.dir/depend
