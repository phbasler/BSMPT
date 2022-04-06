# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later


function(export_headers TARGET HEADER_SOURCE_DIR HEADER_DEST_DIR)

    # Put all headers that are in the source directory into EXPORT_HEADERS variable
    file(GLOB_RECURSE EXPORT_HEADERS CONFIGURE_DEPENDS
        RELATIVE "${HEADER_SOURCE_DIR}"
        "${HEADER_SOURCE_DIR}/*.h"
    )

    # For each header that will be exported
    foreach(HEADER ${EXPORT_HEADERS})

        # Get the directory portion that needs to be created
        get_filename_component(HEADER_DIRECTORY "${HEADER}" DIRECTORY)

        # Create the directory
        add_custom_command(TARGET ${TARGET} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory "${HEADER_DEST_DIR}/${HEADER_DIRECTORY}"
        )

        if (MSVC)

            # Make a hard link to the file
            add_custom_command(TARGET ${TARGET} POST_BUILD
                COMMAND if not exist "${HEADER_DEST_DIR}/${HEADER}" \( mklink /h "${HEADER_DEST_DIR}/${HEADER}" "${HEADER_SOURCE_DIR}/${HEADER}" \)
            )

        else()

            # Make a symbolic link to the file
            add_custom_command(TARGET ${TARGET} POST_BUILD
                COMMAND ln -sf "${HEADER_SOURCE_DIR}/${HEADER}" "${HEADER_DEST_DIR}/${HEADER}"
            )

        endif()


    endforeach(HEADER)

endfunction()


set(header
  C2HDM.h
  R2HDM.h
  RN2HDM.h
  CPINTHEDARK.h
  CXSM.h
  )

set(src
  C2HDM.cpp
  R2HDM.cpp
  RN2HDM.cpp
  CPINTHEDARK.cpp
  CXSM.cpp
  )

add_library(TestCompares ${header} ${src})
target_compile_features(TestCompares PUBLIC cxx_std_14)
export_headers(TestCompares ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/externalincludes/GenerateTestCompares)
target_link_libraries(TestCompares Minimizer)
target_include_directories(TestCompares INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/externalincludes)



