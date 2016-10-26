# 環境変数からboost, python, eigenのインストール先を調べるので、例えば以下のようなコマンドを打つこと。
# export BOOST_ROOT=/opt/boost/
# export NUMPY_INCLUDE_DIR=/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include/
# cmake .

# CMakeのバージョン
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# プロジェクト名
SET(PROJNAME "icassp2010")
PROJECT(${PROJNAME})
# バージョン番号
SET(serial "0.0.1")

# python
#find_package(PythonLibs REQUIRED) # for ${PYTHON_LIBRARY} ${PYTHON_INCLUDE_DIRS}
#MESSAGE(${PYTHON_INCLUDE_DIRS})
#IF(NOT EXISTS ${PYTHON_INCLUDE_DIRS}/numpy/arrayobject.h)
#    IF(NOT DEFINED ENV{NUMPY_INCLUDE_DIR} )
#        MESSAGE( FATAL_ERROR "Please set environment variable NUMPY_INCLUDE_DIR to the include directory of your numpy/arrayobject.h installation.")
#    ENDIF()
#ENDIF()
SET(PYTHON_LIBRARY      "~/.pyenv/versions/anaconda2-4.1.1/lib/")
SET(PYTHON_INCLUDE_DIRS "~/.pyenv/versions/anaconda2-4.1.1/include/python2.7/")
SET(NUMPY_INCLUDE_DIR   "~/.pyenv/versions/anaconda2-4.1.1/lib/python2.7/site-packages/numpy/core/include/")


# build options (currently only for mac)
SET(CMAKE_CXX_FLAGS_DEBUG   "--std=c++11 --stdlib=libc++ -Wall -Wno-unused-function -g ")
SET(CMAKE_CXX_FLAGS_RELEASE "--std=c++11 --stdlib=libc++ -Wall -Wno-unused-function -O2")
# ビルドのタイプはリリース
SET(CMAKE_BUILD_TYPE Release)

# ビルドするためにインクルードするヘッダーファイルの場所(PROJECT_SOURCE_DIRはCMakeList.txtのディレクトリ)
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}" ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR})

# ライブラリの検索パスを指定(複数指定可)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_LIBRARY})

# 実行ファイルを作るなら add_executable(実行ファイル名 リンクさせるファイル...)
# スタティックライブラリを作るなら add_library(ライブラリ名 　STATIC 　リンクさせるファイル...)
# 共有ライブラリ(so)にするなら add_library(ライブラリ名 　SHARED 　リンクさせるファイル...)
FUNCTION(ADD_CYTHON pyxsource)
    # change ext from pyx to cpp
    GET_FILENAME_COMPONENT(name_we ${pyxsource} NAME_WE  )
    GET_FILENAME_COMPONENT(dir     ${pyxsource} DIRECTORY)
    SET(cppsource "${dir}/${name_we}.cpp")
    
    ADD_LIBRARY          (${name_we} MODULE ${cppsource})
    TARGET_LINK_LIBRARIES(${name_we} -lpython2.7)
    ADD_CUSTOM_COMMAND   (OUTPUT ${cppsource}  COMMAND cython ${pyxsource} --cplus  DEPENDS ${pyxsource})
    SET_TARGET_PROPERTIES(${name_we} PROPERTIES PREFIX "" SUFFIX .so LIBRARY_OUTPUT_DIRECTORY_RELEASE ${dir})
ENDFUNCTION()

ADD_CYTHON(prg/signalproc/higmm_core.pyx)