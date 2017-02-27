#ifndef WKY_SIMPLIFIER_H
#define WKY_SIMPLIFIER_H

#include <string>
#include <zjucad/matrix/matrix.h>

namespace SBV {
    using matrixr_t = zjucad::matrix<double>;
    using matrixs_t = zjucad::matrix<size_t>;

    struct Mesh{
        matrixr_t vertices;
        matrixs_t triangles;
    };

    class Simplifier{
    public:
        Simplifier(const std::string& inputMeshPath, const std::string& outputDirectory, double errorBound);

    private:
        std::string mOutputDirectory;
        Mesh mSourceMesh;
    };
}

#endif
