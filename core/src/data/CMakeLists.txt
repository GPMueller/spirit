set(SOURCE_SPIRIT_DATA
    ${SOURCE_SPIRIT_DATA}
    ${CMAKE_CURRENT_SOURCE_DIR}/Geometry.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Spin_System.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Spin_System_Chain.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Spin_System_Chain_Collection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Parameters_Method.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Parameters_Method_LLG.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Parameters_Method_GNEB.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Parameters_Method_MC.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Parameters_Method_MMF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    PARENT_SCOPE # needed so the change of ${SOURCE} will persist to the parent scope
)