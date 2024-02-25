MACRO(sosim_install_all LIB)
    install(TARGETS ${LIB}
            RUNTIME DESTINATION ${SOSIM_INSTALL_DIR}/bin
            LIBRARY DESTINATION ${SOSIM_INSTALL_DIR}/lib
            ARCHIVE DESTINATION ${SOSIM_INSTALL_DIR}/lib)
ENDMACRO()

MACRO(sosim_install_dependency LIBS TARGET)
    if (MSVC)
        add_custom_command(TARGET ${LIB}
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:${LIB}>
                $<TARGET_FILE_DIR:${TARGET}>)
    endif ()
ENDMACRO()

MACRO(sosim_set_exe_output TARGET)
    if (MSVC)
        set_target_properties(${TARGET} PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${SOSIM_INSTALL_DIR}/bin
                LIBRARY_OUTPUT_DIRECTORY ${SOSIM_INSTALL_DIR}/bin
                ARCHIVE_OUTPUT_DIRECTORY ${SOSIM_INSTALL_DIR}/lib
        )
    endif ()
ENDMACRO()