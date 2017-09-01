#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <wkylib/mesh/IO.h>
#include <wkylib/CmdLine.h>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <limits>
#include <eigen3/Eigen/Eigen>
#include <map>
#include "matrix_mn.h"

using namespace zjucad::matrix;
using namespace WKYLIB::Mesh;

matrixr_t points;
matrixs_t lines;

typedef zjucad::matrix::matrix_mn<double, 2, 1> vec2_t;
typedef zjucad::matrix::matrix_mn<double, 3, 1> vec3_t;
using Color = vec3_t;

enum GenType{
    GEN_RESULT = 1,
    GEN_G = 2,
    GEN_GN = 3
};

struct Camera{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

const int g_sampleCount = 500;
const double PI = acos(-1.0);
const double EPSILON = 1e-6;
const int xRes = 400;
const int yRes = 400;
double g_color[xRes][yRes] = {0};
double g_sampleLength;
double cl = -1;
double cr = 1;
matrixr_t un;
std::vector<vec2_t> g_samples;
std::vector<vec2_t> g_samplesN;
int g_totalCount = 0;
GenType g_genType = GEN_RESULT;

inline double distance(const vec2_t& x, const vec2_t& x2)
{
    const double r = norm(x - x2);
    if(fabs(r) < 1e-6) {
        std::cerr << "# warning: near the boundary " << r << std::endl;
    }
    return r;
}

inline double Gn(const vec2_t& x, const vec2_t& x2, const vec2_t& n)
{
    const vec2_t xx = x - x2;
    return (xx[0] * n[0] + xx[1] * n[1]) / (2 * PI * distance(x, x2) * distance(x, x2));
}

inline double G(const vec2_t& x, const vec2_t& x2)
{
    return log(distance(x, x2)) / (2 * PI);
}

inline double E(int sampleIndex)
{
    return un[sampleIndex] + un[sampleIndex + g_totalCount];
}

double computeIntegral(const vec2_t& x, const vec2_t& x2, const vec2_t& n, int sampleIndex)
{
    double result;
    switch(g_genType)
    {
    case GEN_RESULT:
        result = (cr - cl) * Gn(x, x2, n) - G(x, x2) * E(sampleIndex);
        break;
    case GEN_G:
        result = -G(x, x2) * E(sampleIndex);
        break;
    case GEN_GN:
        result = (cr - cl) * Gn(x, x2, n);
        break;
    default:
        throw std::runtime_error("invalid genType.");
    }
    return result;
}

double computeColor(const vec2_t& pixel)
{
    double result = 0;
    for(int i = 0; i < g_totalCount; i++)
    {
        result += computeIntegral(pixel, g_samples[i], g_samplesN[i], i) * g_sampleLength;
    }
    return result;
}

double computeLength()
{
    double result = 0;
    for(int i = 0; i < lines.size(2); i++)
    {
        int a = lines(0, i);
        int b = lines(1, i);
        result += norm(points(colon(), a) - points(colon(), b));
    }
    return result;
}

void computeUN()
{
    matrixr_t A = zeros<double>(g_totalCount * 2, g_totalCount * 2);
    for(int i = 0; i < g_totalCount; i++)
    {
        for(int j = 0; j < g_totalCount; j++)
        {
            if(j == i)
                continue;

            double temp = -G(g_samples[i], g_samples[j]);
            A(i, j) = temp;
            A(i, j + g_totalCount) = temp;
            A(i + g_totalCount, j) = temp;
            A(i + g_totalCount, j + g_totalCount) = temp;
        }
    }
    matrixr_t B = zeros<double>(g_totalCount * 2, 1);
    for(int i = 0; i < g_totalCount; i++)
    {
        B[i] += cl;
        B[i + g_totalCount] += cr;
        for(int j = 0; j < g_totalCount; j++)
        {
            if(j == i)
                continue;

            double temp = cl * Gn(g_samples[i], g_samples[j], -g_samplesN[j]) + cr * Gn(g_samples[i], g_samples[j], g_samplesN[j]);
            B[i] += temp;
            B[i + g_totalCount] += temp;
        }
    }
    std::cout << "solving un..." << std::endl;
 //   Eigen::MatrixXd AA = Eigen::Map<Eigen::MatrixXd>(&A[0], A.size(1), A.size(2));
 //   Eigen::VectorXd bb = Eigen::Map<Eigen::VectorXd>(&B[0], B.size());
 //   Eigen::FullPivLU<Eigen::MatrixXd> solver;
 //   solver.compute(AA);
 //   Eigen::VectorXd dx = solver.solve(bb);
    un = inv(A) * B;
    std::cout << "un solved." << std::endl;
}

void computeAABB(double& xmin, double& xmax, double& ymin, double& ymax)
{
    for(int i = 0; i < points.size(2); i++)
    {
        vec2_t p = points(colon(), i);
        xmin = std::min(xmin, p[0]);
        xmax = std::max(xmax, p[0]);
        ymin = std::min(ymin, p[1]);
        ymax = std::max(ymax, p[1]);
    }
}

void saveImage()
{
    cv::Mat image(xRes, yRes, CV_8UC3, cv::Scalar(0, 0, 0));
    for(int i = 0; i < xRes; i++)
    {
        for(int j = 0; j < yRes; j++)
        {
            cv::Vec3b c;
            if(g_color[i][j] > 0)
            {
                c.val[2] = g_color[i][j] * 255;
                c.val[1] = c.val[0] = 0;
            }
            else{
                c.val[0] = -g_color[i][j] * 255;
                c.val[1] = c.val[2] = 0;
            }
            //const double steps[] = {-0.5, -0.3, -0.2, -0.1, -0.03, 0, 0.03, 0.1, 0.2, 0.3, 0.5};
            const double steps[] = {0.1, 0.3, 0.5, 0.7, 0.9};
            for(size_t k = 0; k < sizeof(steps)/sizeof(double); ++k) {
                if(fabs(fabs(g_color[i][j]) - fabs(steps[k])) < 5e-3)
                {
                    c.val[1] = 255;
                    c.val[0] = c.val[2] = 0;
                }
            }
            image.at<cv::Vec3b>(cv::Point(i, j)) = c;
        }
    }
    std::vector<int> compression_params;
    compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    cv::imwrite("output.png", image, compression_params);
}

void saveGrayImage()
{
    cv::Mat image(xRes, yRes, CV_8UC3, cv::Scalar(0, 0, 0));
    for(int i = 0; i < xRes; i++)
    {
        for(int j = 0; j < yRes; j++)
        {
            cv::Vec3b c;
            c.val[0] = c.val[1] = c.val[2] = g_color[i][j] > 0 ? g_color[i][j] * 255 : 0;
            image.at<cv::Vec3b>(cv::Point(i, j)) = c;
        }
    }
    std::vector<int> compression_params;
    compression_params.push_back(CV_IMWRITE_JPEG_QUALITY);
    cv::imwrite("output.jpg", image, compression_params);
}

void buildAdjacency(std::vector<std::vector<int>>& adjacency)
{
    adjacency.resize(points.size(2));
    for(int i = 0; i < lines.size(2); i++)
    {
        int a = lines(0, i), b = lines(1, i);
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
    }
}

void adaptCurve()
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

void generateSamples()
{
    adaptCurve();

    double currentLength = 0;
    for(int i = 0; i < lines.size(2); i++)
    {
        vec2_t a = points(colon(), lines(0, i));
        vec2_t b = points(colon(), lines(1, i));
        double length = norm(a - b);
        vec2_t nAB = (b - a) / length;
        vec2_t n;
        n[0] = nAB[1];
        n[1] = -nAB[0];
        while(currentLength < length - EPSILON)
        {
            g_samples.push_back(a + currentLength * nAB);
            g_samplesN.push_back(n);
            g_totalCount++;
            currentLength += g_sampleLength;
        }
        currentLength -= length;
    }
    std::cout << "total sample count: " << g_totalCount << std::endl;
    matrixr_t sam(2, g_totalCount);
    for(int i = 0; i < g_totalCount; i++)
    {
        sam(colon(), i) = g_samples[i];
    }
    WKYLIB::Mesh::writePoints2D("samples.vtk", sam);
}

Camera generateCamera()
{
    const double scale = 1.2;

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();

    computeAABB(xmin, xmax, ymin, ymax);

    double xCenter = (xmin + xmax) / 2;
    double yCenter = (ymin + ymax) / 2;
    double xLength = xmax - xmin;
    double yLength = ymax - ymin;
    xLength *= scale;
    yLength *= scale;

    Camera result;
    result.xmin = xCenter - xLength / 2;
    result.xmax = xCenter + xLength / 2;
    result.ymin = yCenter - yLength / 2;
    result.ymax = yCenter + yLength / 2;
    return result;
}

void generateImage()
{
    generateSamples();   
    computeUN();

    Camera camera = generateCamera();
    double xInc = (camera.xmax - camera.xmin) / xRes;
    double yInc = (camera.ymax - camera.ymin) / yRes;

#pragma omp parallel for
    for(int i = 0; i < xRes; i++)
    {
        if(i%10 == 0)
            std::cout << "computing line " << i << "." << std::endl;
        for(int j = 0; j < yRes; j++)
        {
            vec2_t x;
            x[0] = camera.xmin + xInc * i;
            x[1] = camera.ymax - yInc * j;
            g_color[i][j] = computeColor(x);
        }
    }
    double max = *std::max_element(&g_color[0][0], &g_color[0][0]+xRes*yRes),
            min = *std::min_element(&g_color[0][0], &g_color[0][0]+xRes*yRes);
    std::cout << "min max: " << min << " " << max << std::endl;

    for(int i = 0; i < xRes; i++)
        for(int j = 0; j < yRes; j++)
            if(fabs(g_color[i][j]) > 1)
                g_color[i][j] = g_color[i][j] > 0 ? 1 : -1;

    saveImage();
}

int main(int argc, char** argv)
{   
    WKYLIB::Mesh::readCurve2D(argv[1], points, lines);
    g_genType = static_cast<GenType>(std::atoi(argv[2]));

    double length = computeLength();
    g_sampleLength = length / g_sampleCount;

    generateImage();

    return 0;
}
