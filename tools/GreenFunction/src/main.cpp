#include "Common.h"
#include <zjucad/matrix/io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <wkylib/mesh/IO.h>
#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <limits>
#include <map>
#include "Curve.h"
#include "SamplePoint.h"

using namespace zjucad::matrix;
using namespace WKYLIB::Mesh;

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
const double FARAWAY = 2;
const int xRes = 400;
const int yRes = 400;
Curve g_curve;
double g_color[xRes][yRes] = {0};
double g_sampleLength;
double cl = -1;
double cr = 1;
matrixr_t un;
std::vector<SamplePoint> g_samples;
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
    //return 1 / (2 * PI) * distance(x, x2);
}

double computeIntegral(const vec2_t& x, const SamplePoint& sample)
{
    double result;
    const vec2_t& x2 = sample.position;
    const vec2_t& n = sample.normal;
    switch(g_genType)
    {
    case GEN_RESULT:
        result = sample.color * Gn(x, x2, n) - G(x, x2) * sample.derivative;
        break;
    case GEN_G:
        result =  -G(x, x2) * sample.derivative;
        break;
    case GEN_GN:
        result = sample.color * Gn(x, x2, n);
        break;
    default:
        throw std::runtime_error("invalid genType.");
    }
    return result;
}

double computeColor(const vec2_t& pixel)
{
    double result = 0;
    for(int i = 0; i < g_samples.size(); i++)
    {
        result += computeIntegral(pixel, g_samples[i]) * g_sampleLength;
    }
    return result;
}

double computeLength()
{
    double result = 0;
    for(int i = 0; i < g_curve.lines.size(2); i++)
    {
        int a = g_curve.lines(0, i);
        int b = g_curve.lines(1, i);
        result += norm(g_curve.points(colon(), a) - g_curve.points(colon(), b));
    }
    return result;
}

bool isOpposite(const SamplePoint& a, const SamplePoint& b)
{
    return (norm(a.position - b.position) < 1e-6) && (norm(a.normal + b.normal) < 1e-6);
}

void computeUN()
{
    const double l = g_sampleLength;
    const size_t N = g_samples.size();

    matrixr_t A = zeros<double>(N, N);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(j == i || isOpposite(g_samples[i], g_samples[j]))
                A(i, j) = -(-l + l * log(l / 2)) / (2 * PI);
            else
                A(i, j) = -G(g_samples[i].position, g_samples[j].position) * l;
        }
    }

    matrixr_t B = zeros<double>(N, 1);
    for(int i = 0; i < N; i++)
    {
        B[i] += g_samples[i].color;
        for(int j = 0; j < N; j++)
        {
            if(j == i)
                B[i] -= 0.5 * g_samples[j].color;
            else if(isOpposite(g_samples[i], g_samples[j]))
                B[i] += 0.5 * g_samples[j].color;
            else
                B[i] -= (g_samples[j].color * Gn(g_samples[i].position, g_samples[j].position, g_samples[j].normal)) * l;
        }
    }

    std::cout << "solving un..." << std::endl;
    Eigen::Map<Eigen::MatrixXd> AA(&A.data()[0], A.size(1), A.size(2));
    Eigen::Map<Eigen::VectorXd> BB(&B.data()[0], B.size(1), B.size(2));
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(AA);
    Eigen::VectorXd UNN = solver.solve(BB);
    un.resize(UNN.rows(), 1);
    for(int i = 0; i < UNN.rows(); i++){
        un[i] = UNN[i];
    }
    std::cout << un << std::endl;
    std::cout << "un solved." << std::endl;
    for(int i = 0; i < N; i++)
    {
        g_samples[i].derivative = un[i];
    }

    auto lefttop = A(colon(0, N / 2 - 1), colon(0, N / 2 - 1));
    auto leftbottom = A(colon(N / 2, N - 1), colon(0, N / 2 - 1));
}

void computeAABB(double& xmin, double& xmax, double& ymin, double& ymax)
{
    for(int i = 0; i < g_curve.points.size(2); i++)
    {
        vec2_t p = g_curve.points(colon(), i);
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
                c.val[2] = g_color[i][j] > 1 ? 255 : g_color[i][j] * 255;
                c.val[1] = c.val[0] = 0;
            }
            else{
                c.val[0] = g_color[i][j] < -1 ? 255 : -g_color[i][j] * 255;
                c.val[1] = c.val[2] = 0;
            }

            const double steps[] = {0.1, 0.3, 0.5, 0.7, 0.9, 0.99};
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

void generateSamples()
{
    g_curve.adjust();

    double currentLength = 0;
    for(int i = 0; i < g_curve.lines.size(2); i++)
    {
        vec2_t a = g_curve.points(colon(), g_curve.lines(0, i));
        vec2_t b = g_curve.points(colon(), g_curve.lines(1, i));
        double length = norm(a - b);
        vec2_t nAB = (b - a) / length;
        vec2_t n;
        n[0] = nAB[1];
        n[1] = -nAB[0];
        while(currentLength < length - EPSILON)
        {
            SamplePoint sample;
            sample.position = a + currentLength * nAB;
            sample.normal = n;
            sample.color = cr;
            g_samples.push_back(sample);
            sample.normal = -n;
            sample.color = cl;
            g_samples.push_back(sample);
            currentLength += g_sampleLength;
        }
        currentLength -= length;
    }
    std::cout << "total sample count: " << g_samples.size() / 2 << std::endl;
    matrixr_t sam(2, g_samples.size() / 2);
    int j = 0;
    for(int i = 0; i < g_samples.size(); i++)
    {
        if(i % 2)
            sam(colon(), j++) = g_samples[i].position;
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

void initFarAway()
{
    const int N = FARAWAY / g_sampleLength;
    for(int i = 0; i < N; i++)
    {
        const double temp = i * g_sampleLength;
        //right top vertical
        SamplePoint point;
        point.position[0] = FARAWAY;
        point.position[1] = temp;
        point.normal[0] = -1;
        point.normal[1] = 0;
        g_samples.push_back(point);

        //right bottom vertical
        if(i)
        {
            point.position[0] = FARAWAY;
            point.position[1] = -temp;
            point.normal[0] = -1;
            point.normal[1] = 0;
            g_samples.push_back(point);
        }

        //left top vertical
        point.position[0] = -FARAWAY;
        point.position[1] = temp;
        point.normal[0] = 1;
        point.normal[1] = 0;
        g_samples.push_back(point);

        if(i)
        {
            //left bottom vertical
            point.position[0] = -FARAWAY;
            point.position[1] = -temp;
            point.normal[0] = 1;
            point.normal[1] = 0;
            g_samples.push_back(point);
        }

        //right top horizontal
        point.position[0] = temp;
        point.position[1] = FARAWAY;
        point.normal[0] = 0;
        point.normal[1] = -1;
        g_samples.push_back(point);

        //right bottom horizontal
        point.position[0] = temp;
        point.position[1] = -FARAWAY;
        point.normal[0] = 0;
        point.normal[1] = 1;
        g_samples.push_back(point);

        if(i)
        {
            //left top horizontal
            point.position[0] = -temp;
            point.position[1] = FARAWAY;
            point.normal[0] = 0;
            point.normal[1] = -1;
            g_samples.push_back(point);

            //left bottom horizontal
            point.position[0] = -temp;
            point.position[1] = -FARAWAY;
            point.normal[0] = 0;
            point.normal[1] = 1;
            g_samples.push_back(point);
        }
    }
}

void addInternalPoint(const Camera& camera)
{
    vec2_t center;
    center[0] = (camera.xmin + camera.xmax) / 2;
    center[1] = (camera.ymin + camera.ymax) / 2;

    const double radius = 0.01;
    double theta = g_sampleLength / radius;
    for(int i = 0; i * theta < 2 * PI; i++)
    {
        double t = i * theta;
        SamplePoint point;
        point.position[0] = center[0] + radius * cos(t);
        point.position[1] = center[1] + radius * sin(t);
        point.normal[0] = cos(t);
        point.normal[1] = sin(t);
        point.color = 0;
        g_samples.push_back(point);
    }
}

void generateImage()
{
    Camera camera = generateCamera();
    double xInc = (camera.xmax - camera.xmin) / xRes;
    double yInc = (camera.ymax - camera.ymin) / yRes;

    generateSamples();
    addInternalPoint(camera);
    initFarAway();
    computeUN();

#ifdef NDEBUG
#pragma omp parallel for
#endif
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

    // for(int i = 0; i < xRes; i++)
    //    for(int j = 0; j < yRes; j++)
            //g_color[i][j] /= fabs(min);
            //if(fabs(g_color[i][j]) > 1)
           //     g_color[i][j] = g_color[i][j] > 0 ? 1 : -1;


    saveImage();
}



int main(int argc, char** argv)
{   
    WKYLIB::Mesh::readCurve2D(argv[1], g_curve.points, g_curve.lines);
    g_genType = static_cast<GenType>(std::atoi(argv[2]));

    double length = computeLength();
    g_sampleLength = length / g_sampleCount;

    generateImage();

    return 0;
}
