file(GLOB subheaders *.h )
file(GLOB subsources *.cpp )

set(headers
    ${headers}
    ${subheaders}
    PARENT_SCOPE
    )

set(sources
    ${sources}
    ${subsources}
    PARENT_SCOPE
    )
