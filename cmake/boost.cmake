set(BOOST_IOSTREAMS_ENABLE_ZSTD OFF CACHE BOOL "" FORCE)

set(BOOST_IOSTREAMS_ENABLE_LZMA OFF)
set(BOOST_IOSTREAMS_ENABLE_BZIP2 OFF)

include(FetchContent)

if(WIN32)
    FetchContent_Declare(
        boost
        URL https://github.com/boostorg/boost/releases/download/boost-1.92.0/boost-1.92.0-cmake.zip
        URL_HASH SHA256=f47ee00d5870f006ef95a7be307be6553a769428d48ee98f20e856b9c3dfd3a2
        OVERRIDE_FIND_PACKAGE
    )
else()
    FetchContent_Declare(
        boost
        URL https://github.com/boostorg/boost/releases/download/boost-1.92.0/boost-1.92.0-cmake.tar.gz
        URL_HASH SHA256=f51707c27359a0df0cac1beada86de31bb5eed5e8285592dadec384df99c2984
        OVERRIDE_FIND_PACKAGE
    )
endif()

FetchContent_MakeAvailable(boost)

if(zlib_SOURCE_DIR)
    target_include_directories(boost_iostreams PRIVATE ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
endif()
