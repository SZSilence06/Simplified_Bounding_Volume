#include "Curve.h"

using namespace zjucad::matrix;

void Curve::buildAdjacency(std::vector<std::vector<int>>& adjacency)
{
    adjacency.resize(points.size(2));
    for(int i = 0; i < lines.size(2); i++)
    {
        int a = lines(0, i), b = lines(1, i);
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
    }
}

void Curve::adjust()
{
    std::vector<std::vector<int>> adjacency;
    buildAdjacency(adjacency);
    int endPoints[2];
    int k = 0;
    for(int i = 0; i < adjacency.size(); i++)
    {
        if(adjacency[i].size() == 1)
        {
            endPoints[k++] = i;
        }
    }
    if(k == 0)
    {
        endPoints[0] = endPoints[1] = 0;
    }
    k = 0;
    int prev = -1;
    for(int i = endPoints[0]; i != endPoints[1];)
    {
        for(int j = 0; j < adjacency[i].size(); j++)
        {
            if(adjacency[i][j] != prev)
            {
                matrixs_t line(2, 1);
                line[0] = i;
                line[1] = adjacency[i][j];
                lines(colon(), k++) = line;
                prev = i;
                i = adjacency[i][j];
            }
        }
    }
}
