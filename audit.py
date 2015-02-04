#!/usr/bin/env python
#
# Copyright 2013 German Aerospace Center (http://www.DLR.de)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#       http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This script checks for common C++ programming isssues and coding standard.
It uses the tools cppcheck and vera++ to perform these checks. It is expected that both
command line tools are available in the PATH environment.

Usage:
* python audit.py
 * Prints the results on standard output
 * Should be run before you commit
* python audit.py xml
 * Captures results in cppcheck.xml and vera++.xml
 * Is intended for the continuous integration build
"""


import fnmatch
import os
import subprocess
import sys


_XML_OUT_OPTION = "xml"
_VERA_BASE_PATH = os.path.abspath(os.path.join(os.curdir, "tools", "style-check"))
_PROJECT_ROOT_PATH = os.path.abspath(os.curdir)

# Names of the top-level directories which should not be analyzed
_EXCLUDED_TOP_LEVEL_DIRECTORIES = [".git", "build", "config", "doc", "examples", "tools", "CMakeFiles", "_CPack_Packages"]

# Additional include path (for cppcheck)
_INCLUDE_PATH = "include"

def audit(outputformat=None):
    """ Controls the check process.

    @param outputformat: Indicates the overall output format.
        xml => Produces XML reports
        everything else => Prints results on standard output
    """

    excluded_paths = _determine_excluded_paths()
    headers = _find_files("*.hpp", excluded_paths)
    sources = _find_files("*.cpp", excluded_paths)

    _run_cppcheck(outputformat, headers, sources)
    _run_vera(outputformat, headers, sources)


def _determine_excluded_paths():
    """ Just builds the absolute paths of the excluded top-level directories. """

    excluded_dirs = list()
    for directory_name in _EXCLUDED_TOP_LEVEL_DIRECTORIES:
        excluded_dirs.append(os.path.join(_PROJECT_ROOT_PATH, directory_name))
    return excluded_dirs


def _find_files(pattern, excluded_paths):
    """ Returns a list of all file paths which match the given
    Unix style filename pattern. It recursively walks all directories
    of the current directory.
    """

    matches = list()
    for current_path, _, filenames in os.walk(_PROJECT_ROOT_PATH):
        if not _exclude(current_path, excluded_paths):
            for filename in fnmatch.filter(filenames, pattern):
                matches.append(os.path.join(current_path, filename))
    return matches

def _exclude(current_path, excluded_paths):
    """ Checks whether the actual path needs to be excluded. """

    for excluded_path in excluded_paths:
        if current_path == excluded_path:
            return True
        if current_path.startswith(excluded_path):
            if current_path[len(excluded_path)] == os.sep:
                return True
    return False


def _run_cppcheck(outputformat, headers, sources):
    """ Runs cppcheck with all available checks and uses the C++11 profile. """

    # Create the command
    cppcheck_command = "cppcheck --std=c++03 --enable=information,performance,portability,style"
    cppcheck_command += " -I " + _INCLUDE_PATH
    if outputformat == _XML_OUT_OPTION:
        cppcheck_command += " --xml --xml-version=2"
    else:
        cppcheck_command += " --template=\"{file}:{line}: warning: cppcheck:{message}\""
    files_tocheck = " ".join(headers)
    files_tocheck += " " + " ".join(sources)
    cppcheck_command += " %s" % files_tocheck

    # Run the command
    _, _, stderr = _run_command(cppcheck_command)

    # Correctly handle the result
    if outputformat == _XML_OUT_OPTION:
        file_object = open("cppcheck.xml", "wb")
        file_object.write(stderr)
        file_object.close()
    else:
        print(stderr)


def _run_vera(outputformat, headers, sources):
    """ Runs vera++ with the self-defined rules. """

    # Create the command
    files_tocheck = " ".join(headers)
    files_tocheck += " " + " ".join(sources)
    vera_command = "vera++ %s" % files_tocheck
    #vera_command += " -r %s" % _VERA_BASE_PATH # Vera root directory to find the right scripts/rules/profiles
    vera_command += " --warning" # Just sho gcc warnings
    vera_command += " -p gm" # Set the right profile
    #vera_command += " --exclusions %s" % os.path.join(_VERA_BASE_PATH, "scripts", "exclusions.tcl")
    if outputformat == _XML_OUT_OPTION:
        vera_command += " -c vera++.xml" # Optional XML checkstyle output for Jenkins

    # Run the command
    _, _, stderr = _run_command(vera_command)

    # Optionally print the putput on standard output
    if outputformat != _XML_OUT_OPTION:
        print(stderr)


def _run_command(command):
    """ Runs the given command and returns exit code, standard out/error. """

    print("Running command: %s" % command)
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    exit_code = process.returncode
    if exit_code != 0:
        print(stdout, stderr)
    return exit_code, stdout, stderr


if __name__ == "__main__":
    outputformat = ""
    if len(sys.argv) == 2:
        outputformat = sys.argv[1]
    audit(outputformat=outputformat)
