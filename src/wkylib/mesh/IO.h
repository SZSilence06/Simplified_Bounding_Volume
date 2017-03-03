#ifndef WKY_MESH_IO_H
#define WKY_MESH_IO_H

#include "Common.h"

namespace WKYLIB
{
    namespace Mesh
    {
        bool writePoints(const matrixr_t& points, const std::string& file);
    }
}

#endif
