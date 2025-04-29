# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

add_compile_options(
  $<$<CONFIG:DEBUG>:-DCOMPILEDEBUG=true>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-pedantic>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wshadow>
  $<$<AND:$<CONFIG:DEBUG>,$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>>:-Wmissing-declarations>
  $<$<AND:$<CONFIG:DEBUG>,$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>>:-Wmissing-include-dirs>
  $<$<AND:$<CONFIG:RELEASE>,$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>>:-O3>
)

if(BSMPTUseVectorization)
  include(CheckCXXCompilerFlag)

  check_cxx_compiler_flag("-march=native" _march_native_works)
  check_cxx_compiler_flag("-xHost" _xhost_works)

  if(_march_native_works)
    message(
      STATUS
        "Using processor's vector instructions (-march=native compiler flag set)"
    )
    add_compile_options(-march=native)
  elseif(_xhost_works)
    message(
      STATUS "Using processor's vector instructions (-xHost compiler flag set)")
    add_compile_options(-xHost)
  else()
    message(STATUS "No suitable compiler flag found for vectorization")
  endif()
endif(BSMPTUseVectorization)

add_compile_options(
  $<$<CXX_COMPILER_ID:MSVC>:/permissive->
  $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:DEBUG>>:/bigobj>
  $<$<CXX_COMPILER_ID:MSVC>:/w44101>
  $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:RELEASE>>:/Ox>)

list(
  APPEND
  MSVC_DISABLED_WARNINGS_LIST
  "C4061" # enumerator 'identifier' in switch of enum 'enumeration' is not
          # explicitly handled by a case label Disable this because it flags
          # even when there is a default.
  "C4068"
  "C4100" # 'exarg' : unreferenced formal parameter
  "C4127" # conditional expression is constant
  "C4200" # nonstandard extension used : zero-sized array in struct/union.
  "C4204" # nonstandard extension used: non-constant aggregate initializer
  "C4221" # nonstandard extension used : 'identifier' : cannot be initialized
          # using address of automatic variable
  "C4242" # 'function' : conversion from 'int' to 'uint8_t', possible loss of
          # data
  "C4244" # 'function' : conversion from 'int' to 'uint8_t', possible loss of
          # data
  "C4245" # 'initializing' : conversion from 'long' to 'unsigned long',
          # signed/unsigned mismatch
  "C4267" # conversion from 'size_t' to 'int', possible loss of data
  "C4355"
  "C4371" # layout of class may have changed from a previous version of the
          # compiler due to better packing of member '...'
  "C4388" # signed/unsigned mismatch
  "C4296" # '>=' : expression is always true
  "C4350" # behavior change: 'std::_Wrap_alloc...'
  "C4365" # '=' : conversion from 'size_t' to 'int', signed/unsigned mismatch
  "C4389" # '!=' : signed/unsigned mismatch
  "C4464" # relative include path contains '..'
  "C4510" # 'argument' : default constructor could not be generated
  "C4530" # C++ exception handler used
  "C4571"
  "C4512" # 'argument' : assignment operator could not be generated
  "C4514" # 'function': unreferenced inline function has been removed
  "C4548" # expression before comma has no effect; expected expression with
          # side-effect" caused by FD_* macros.
  "C4610" # struct 'argument' can never be instantiated - user defined
          # constructor required.
  "C4619"
  "C4623" # default constructor was implicitly defined as deleted
  "C4625" # copy constructor could not be generated because a base class copy
          # constructor is inaccessible or deleted
  "C4626" # assignment operator could not be generated because a base class
          # assignment operator is inaccessible or deleted
  "C4643"
  "C4668" # 'symbol' is not defined as a preprocessor macro, replacing with '0'
          # for 'directives' Disable this because GTest uses it everywhere.
  "C4706" # assignment within conditional expression
  "C4710" # 'function': function not inlined
  "C4711" # function 'function' selected for inline expansion
  "C4800" # 'int' : forcing value to bool 'true' or 'false' (performance
          # warning)
  "C4820" # 'bytes' bytes padding added after construct 'member_name'
  "C4868"
  "C4996"
  "C5026" # move constructor was implicitly defined as deleted
  "C5027" # move assignment operator was implicitly defined as deleted
  "C5031"
  "C5039"
  "C5045")

foreach(warning IN LISTS MSVC_DISABLED_WARNINGS_LIST)
  string(REPLACE "C" "" warning_stripped ${warning})
  add_compile_options($<$<CXX_COMPILER_ID:MSVC>:-wd${warning_stripped}>)
endforeach()
