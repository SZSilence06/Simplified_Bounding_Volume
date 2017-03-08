#include <wkylib/mesh/IO.h>

using matrixr_t = WKYLIB::Mesh::matrixr_t;
using matrixs_t = WKYLIB::Mesh::matrixs_t;

using namespace zjucad::matrix;

int main(int argc, char** argv)
{
    matrixr_t vertices;
    matrixs_t lines;

    if(WKYLIB::Mesh::readCurve2D(argv[1], vertices, lines) == false)
    {
        std::cout << "fail to read curve file : " + std::string(argv[1]) << std::endl;
        return 0;
    }

    std::vector<int> prior;
    std::vector<int> next;

    prior.reserve(vertices.size(2));
    next.reserve(vertices.size(2));
    for(int i = 0; i < vertices.size(2); i++)
    {
        prior.push_back(-1);
    }
    for(int i = 0; i < vertices.size(2); i++)
    {
        next.push_back(-1);
    }

    for(int i = 0; i < lines.size(2); i++)
    {
        if(prior[lines(1, i)] < 0)
        {
            prior[lines(1, i)] = lines(0, i);
        }
        else
        {
            next[lines(1, i)] = lines(0, i);
        }
        if(next[lines(0, i)] < 0)
        {
            next[lines(0, i)] = lines(1, i);
        }
        else
        {
            prior[lines(0, i)] = lines(1, i);
        }
    }

    int start = lines(0, 0);
    int p1 = start;
    int p2 = next[p1] == -1 ? prior[p1] : next[p1];

    matrixs_t new_lines(2, lines.size(2));
    int i = 0;
    while(p2 != start && i < vertices.size(2))
    {
        new_lines(0, i) = p1;
        new_lines(1, i) = p2;

        if(prior[p2] != p1)
        {
            int temp = prior[p2];
            next[p1] = p2;
            prior[p2] = p1;
            p1 = p2;
            p2 = temp;
        }
        else
        {
            next[p1] = p2;
            p1 = p2;
            p2 = next[p2];
        }

        i++;
    }
    new_lines(0, i) = p1;
    new_lines(1, i) = p2;

    WKYLIB::Mesh::writeCurve2D("adjusted_" + std::string(argv[1]), vertices, new_lines);

    return 0;
}
