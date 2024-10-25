# CMake 최소 요구 버전
cmake_minimum_required(VERSION 3.7)

# 프로젝트 이름 및 C++ 표준 설정
project(pRNE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 툴체인 경로 및 설정을 위한 옵션 플래그
option(BUILD_FOR_COMPACTRIO "Build for CompactRIO (Linux-based)" ON)
option(BUILD_FOR_WINDOWS "Build for Windows" OFF)
option(BUILD_FOR_UBUNTU "Build for Ubuntu" OFF)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
set(LIBNAME ${CMAKE_PROJECT_NAME})
set(LIBNAME_UBUNTU ${CMAKE_PROJECT_NAME}_ubuntu)
if (CMAKE_BUILD_TYPE STREQUAL  "")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")
endif()
# 소스 파일 추가
set(SOURCES src/gen.c)
set(TEST_SOURCES src/callingSharedObject.cpp)
set(TEST_SOURCES_WIN src/callingSharedObject_win.cpp)

# CompactRIO 타겟 설정 (Linux 환경)
if(BUILD_FOR_COMPACTRIO)
    message(STATUS "Building for CompactRIO...")

    # 타겟 시스템 정보 설정 (NI Linux Real-Time x64)
    set(CMAKE_SYSTEM_NAME Linux)
    set(CMAKE_SYSTEM_PROCESSOR x86_64)

    # 툴체인 경로 설정 (Linux 환경)
    set(toolchainpath /usr/local/oecore-x86_64/sysroots)

    # 컴파일러 설정
    set(CMAKE_C_COMPILER ${toolchainpath}/x86_64-nilrtsdk-linux/usr/bin/x86_64-nilrt-linux/x86_64-nilrt-linux-gcc)
    set(CMAKE_CXX_COMPILER ${toolchainpath}/x86_64-nilrtsdk-linux/usr/bin/x86_64-nilrt-linux/x86_64-nilrt-linux-g++)

    # Sysroot 설정
    set(CMAKE_SYSROOT ${toolchainpath}/core2-64-nilrt-linux)

    # 컴파일러 플래그 설정
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fmessage-length=0")
    
    # Release 모드 설정
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")  # Release 모드 최적화 플래그
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")  # Debug 모드 플래그
    endif()

    # 인클루드 디렉토리 설정 (헤더 파일 경로)
    include_directories(${CMAKE_SYSROOT}/usr/include)
    include_directories(${CMAKE_SOURCE_DIR}/src)  # src 디렉토리 추가
    include_directories(${CMAKE_SYSROOT}/usr/include/c++/6.3.0)
    include_directories(${CMAKE_SYSROOT}/usr/include)
    include_directories(${CMAKE_SYSROOT}/usr/include/c++/6.3.0/x86_64-nilrt-linux)
    
    # 링크 디렉토리 설정
    link_directories(${CMAKE_SYSROOT}/usr/lib)

    # RPATH 설정 (런타임에 라이브러리 경로)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_SYSROOT}/usr/lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # CompactRIO용 `.so` 라이브러리 생성
    add_library(${LIBNAME} SHARED ${SOURCES})
    target_include_directories(${LIBNAME} PUBLIC 
        ${CMAKE_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
    )
    
    # 전처리기 플래그 추가 (Linux용)
    add_definitions(-DCOMPACTRIO)

    # CompactRIO 빌드 시 다른 타겟 설정 무시
    set(BUILD_FOR_WINDOWS OFF)
    set(BUILD_FOR_UBUNTU OFF)

endif()

# Windows 타겟 설정 (DLL - 32bit MinGW 사용)
if(BUILD_FOR_WINDOWS)
    message(STATUS "Building for Windows (32-bit)...")

    # 32비트 MinGW 컴파일러 설정 (Windows용 크로스 컴파일)
    set(CMAKE_SYSTEM_NAME Windows)
    set(CMAKE_SYSTEM_PROCESSOR i686)  # 32비트 프로세서로 설정
    set(CMAKE_C_COMPILER i686-w64-mingw32-gcc)
    set(CMAKE_CXX_COMPILER i686-w64-mingw32-g++)
    set(CMAKE_SHARED_LINKER_FLAGS "")

    # Windows에서 DLL 파일을 만들기 위한 전처리기 정의
    # 32비트 DLL 파일 생성 (링커에 32비트 타겟 명시)
    set(CMAKE_C_FLAGS "-m32")
    set(CMAKE_CXX_FLAGS "-m32")
    set(CMAKE_EXE_LINKER_FLAGS "-m32")

    # Release 모드 설정
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")  # Release 모드 최적화 플래그
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")  # Debug 모드 플래그
    endif()

    # Windows용 `.dll` 라이브러리 생성
    add_library(${LIBNAME} SHARED ${SOURCES})
    target_include_directories(${LIBNAME} PUBLIC ${CMAKE_SOURCE_DIR}/src)
    include_directories(${CMAKE_SOURCE_DIR}/src)
    include_directories(${PROJECT_SOURCE_DIR}/include)
    
    # DLL 파일 생성 시 필요한 추가 플래그 설정 (MinGW 링커 옵션)
    set_target_properties(${LIBNAME} PROPERTIES PREFIX "" SUFFIX ".dll")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")

    # Windows 테스트 실행 파일 생성
    add_executable(test_win ${TEST_SOURCES_WIN})
    target_link_libraries(test_win ${LIBNAME})
    add_definitions(-DWINDOWS)

    # Windows 빌드 시 다른 타겟 설정 무시
    set(BUILD_FOR_COMPACTRIO OFF)
    set(BUILD_FOR_UBUNTU OFF)
endif()


# Ubuntu 타겟 설정 (일반적인 빌드 방식)
if(BUILD_FOR_UBUNTU)
    message(STATUS "Building for Ubuntu...")

    # 기본 GCC 및 G++ 컴파일러 사용
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)

    # 인클루드 디렉토리 설정
    include_directories(${CMAKE_SOURCE_DIR}/src)
    include_directories(${PROJECT_SOURCE_DIR}/include)

    # Ubuntu용 `.so` 라이브러리 생성
    add_library(${LIBNAME_UBUNTU} SHARED ${SOURCES})
    target_link_libraries(${LIBNAME_UBUNTU} m)

    # 전처리기 플래그 추가 (Ubuntu용)
    add_definitions(-DUBUNTU)
    add_executable(test ${TEST_SOURCES})
    target_link_libraries(test PUBLIC  ${LIBNAME_UBUNTU} dl m)
    # Ubuntu 빌드 시 다른 타겟 설정 무시
    set(BUILD_FOR_COMPACTRIO OFF)
    set(BUILD_FOR_WINDOWS OFF)


endif()

# 타겟이 지정되지 않았을 때 오류 메시지 표시
if(NOT BUILD_FOR_COMPACTRIO AND NOT BUILD_FOR_WINDOWS AND NOT BUILD_FOR_UBUNTU)
    message(FATAL_ERROR "No target specified! Please set either BUILD_FOR_COMPACTRIO, BUILD_FOR_WINDOWS, or BUILD_FOR_UBUNTU to ON.")
endif()
