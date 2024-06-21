#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Open3D::Open3D" for configuration "Release"
set_property(TARGET Open3D::Open3D APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Open3D::Open3D PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/Open3D.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/Open3D.dll"
  )

list(APPEND _cmake_import_check_targets Open3D::Open3D )
list(APPEND _cmake_import_check_files_for_Open3D::Open3D "${_IMPORT_PREFIX}/lib/Open3D.lib" "${_IMPORT_PREFIX}/bin/Open3D.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
