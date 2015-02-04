IF(NOT EXISTS "/home_local/matu_gr/work/cresta/ppstee/trunk/install_manifest.txt")
  MESSAGE(FATAL_ERROR "Cannot find install manifest: "/home_local/matu_gr/work/cresta/ppstee/trunk/install_manifest.txt"")
ENDIF(NOT EXISTS "/home_local/matu_gr/work/cresta/ppstee/trunk/install_manifest.txt")

FILE(READ "/home_local/matu_gr/work/cresta/ppstee/trunk/install_manifest.txt" files)
STRING(REGEX REPLACE "\n" " " files "${files}")
separate_arguments (files)
FOREACH(file ${files})
  MESSAGE(STATUS "Uninstalling "$ENV{DESTDIR}${file})
  IF(EXISTS "$ENV{DESTDIR}${file}")
    EXECUTE_PROCESS(
      COMMAND /tools/modulesystem/tools/cmake/cmake-2.8.11/install/sled11.x86_64.gcc-4.3.4.release/bin/cmake -E remove $ENV{DESTDIR}${file}
      OUTPUT_VARIABLE rm_out
      RESULT_VARIABLE rm_retval
      )
    IF(NOT "${rm_retval}" STREQUAL 0)
      MESSAGE(FATAL_ERROR "Problem when removing "$ENV{DESTDIR}${file})
    ENDIF(NOT "${rm_retval}" STREQUAL 0)
  ELSE(EXISTS "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "File "$ENV{DESTDIR}${file}" does not exist.")
  ENDIF(EXISTS "$ENV{DESTDIR}${file}")
ENDFOREACH(file)
