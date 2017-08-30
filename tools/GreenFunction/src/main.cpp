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

const int g_sampleCount = 5000;
const double PI = acos(-1.0);
const double EPSILON = 5e-3;
const int xRes = 400;
const int yRes = 400;
double g_color[xRes][yRes] = {0};
double g_sampleLength;
double cl = 1;
double cr = 1;
matrixr_t un;
vec2_t g_samples[g_sampleCount];
vec2_t g_samplesN[g_sampleCount];
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
    return (xx[0] * n[0] + xx[1] * n[1]) / (2 * PI * distance(x, x2));
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
    const double scale = 1;
    for(int i = 0; i < g_totalCount; i++)
    {
        result += computeIntegral(pixel*scale, g_samples[i]*scale, g_samplesN[i], i) * g_sampleLength*scale;
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
            const double steps[] = {0.1, 0.5, 0.9};
            for(size_t k = 0; k < sizeof(steps)/sizeof(double); ++k) {
                if(fabs(fabs(g_color[i][j]) - fabs(steps[k])) < EPSILON)
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

void generateSamples()
{
    for(int i = 0; i < lines.size(2); i++)
    {
        vec2_t a = points(colon(), lines(0, i));
        vec2_t b = points(colon(), lines(1, i));
        double length = norm(a - b);
        vec2_t nAB = (b - a) / length;
        vec2_t n;
        n[0] = nAB[1];
        n[1] = -nAB[0];
        int sampleCount = length / g_sampleLength;
        for(int j = 0; j < sampleCount; j++)
        {
            g_samples[g_totalCount] = a + j * g_sampleLength * nAB;
            g_samplesN[g_totalCount] = n;
            g_totalCount++;
        }
    }
    std::cout << "total sample count: " << g_totalCount << std::endl;
}

Camera generateCamera()
{
    const double scale = 1e1;

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
            g_color[i][j] /= max;

    if(g_genType == GEN_RESULT)
    {
        saveGrayImage();
    }
    else
    {
        saveImage();
    }
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
