#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>
#include <GL/glew.h>
#include <GL/freeglut.h>

#define PI 3.14159265

struct Point {
    float x;
    float y;
    float z;
};

typedef std::pair<std::vector<std::string>, std::vector<float>> Transformation;
std::vector<Point> oPoints;
std::vector<Point> pPrimes;

void process();
arma::mat getCompositeMatrix(std::vector<Transformation> &transformations);
std::vector<Point> transform(std::vector<Point> &points, arma::mat &T);
void drawScene();
void resize(int w, int h);
void setup();
Point translate(Point p, float D[]);
Point scale(Point p, float S[]);
Point rotateOnX(Point p, float theta);
Point rotateOnY(Point p, float theta);
Point rotateOnZ(Point p, float theta);
float degToRad(float deg);
void printPoint(Point p);

int main(int argc, char* argv[]) {
    process();

    /* OpenGL related calls. */
    glutInit(&argc, argv);
    glutInitContextVersion(2, 1);
    glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);

    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("transformaciones3d.cpp");

    glutDisplayFunc(drawScene);
    glutReshapeFunc(resize);

    glewInit();

    setup();

    glutMainLoop();

    return EXIT_SUCCESS;
}

void drawScene() {
    glClear(GL_COLOR_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(15.0, 10.0, 15.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    /* Draw red lines to depict the axes. */
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(100.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 100.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 100.0);
    glEnd();

    /* Draw original points. */
    glColor3f(0.0, 0.0, 0.0);
    glBegin(GL_POLYGON);
        for (auto p : oPoints) {
            glVertex3f(p.x, p.y, p.z);
        }
    glEnd();

    /* Draw transformed points. */
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_POLYGON);
        for (auto p : pPrimes) {
            glVertex3f(p.x, p.y, p.z);
        }
    glEnd();

    glFlush();
}

void resize(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-5.0, 5.0, -5.0, 5.0, 10.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void setup() {
    glClearColor(1.0, 1.0, 1.0, 0.0);
}

void process() {
    int n, t;
    arma::mat T;
    while (std::cin >> n) {
        std::cin >> t;
        std::vector<Transformation> transformations;

        /* Read all the points */
        for (int i = 0; i < n; i++) {
            Point p;
            std::cin >> p.x >> p.y >> p.z;
            oPoints.push_back(p);
        }

        /* Read all the vectors of transformations. */
        for (int i = 0; i < t; i++) {
            Transformation t;
            std::string tName;
            std::cin >> tName;
            t.first.push_back(tName);

            /* Treat translation and scaling differently because of the
             * number of parameters */
            if (t.first[0] == "t" || t.first[0] == "s") {
                float c;
                std::cin >> c;
                t.second.push_back(c);
                std::cin >> c;
                t.second.push_back(c);
                std::cin >> c;
                t.second.push_back(c);
            } else if (t.first[0] == "r") {
                float c;
                std::string axis;
                std::cin >> axis;
                t.first.push_back(axis);
                std::cin >> c;
                t.second.push_back(c);
            }
            transformations.push_back(t);
        }
        T = getCompositeMatrix(transformations);
        T.print();
        for (auto r : transform(oPoints, T)) {
            pPrimes.push_back(r);
        }
    }
}

/*
 * Calculates the composite matrix from the vector of individual transformations
 * by multiplying them from left to right.
 */
arma::mat getCompositeMatrix(std::vector<Transformation> &transformations) {
    arma::mat T(4, 4, arma::fill::eye), B(4, 4);

    /* Iterate through all the transformation vectors mapping each to the
     * corresponding transformation matrix in order to perform the product. */
    for (auto t : transformations) {
        if (t.first[0] == "t") {
            B = { {1, 0, 0, t.second[0]},
                  {0, 1, 0, t.second[1]},
                  {0, 0, 1, t.second[2]},
                  {0, 0, 0, 1} };
        } else if (t.first[0] == "s") {
            B = { {t.second[0], 0, 0, 0},
                  {0, t.second[1], 0, 0},
                  {0, 0, t.second[2], 0},
                  {0, 0, 0, 1} };
        } else if (t.first[0] == "r") {
            float theta = degToRad(t.second[0]);
            float cosTheta = cos(theta), sinTheta = sin(theta);
            if (t.first[1] == "x") {
                B = { {1, 0, 0, 0},
                      {0, cosTheta, -sinTheta, 0},
                      {0, sinTheta, cosTheta, 0},
                      {0, 0, 0, 1} };
            } else if (t.first[1] == "y") {
                B = { {cosTheta, 0, sinTheta, 0},
                      {0, 1, 0, 0},
                      {-sinTheta, 0, cosTheta, 0},
                      {0, 0, 0, 1} };
            } else {
                B = { {cosTheta, -sinTheta, 0, 0},
                      {sinTheta, cosTheta, 0, 0},
                      {0, 0, 1, 0},
                      {0, 0, 0, 1} };
            }
        }
        T *= B;
    }
    return T;
}

/*
 * Applies the transformations to each point in the vector.
 */
std::vector<Point> transform(std::vector<Point> &points, arma::mat &T) {
    std::vector<Point> res;
    for (auto point : points) {
        Point p;
        arma::mat v({point.x, point.y, point.z, 1}), pPrime;
        pPrime = T * v.t();
        p.x = pPrime[0];
        p.y = pPrime[1];
        p.z = pPrime[2];
        res.push_back(p);
    }
    return res;
}

Point translate(Point p, float D[]) {
    Point pPrime;
    pPrime.x = p.x + D[0];
    pPrime.y = p.y + D[1];
    pPrime.z = p.z + D[2];
    return pPrime;
}

Point scale(Point p, float S[]) {
    Point pPrime;
    pPrime.x = p.x * S[0];
    pPrime.y = p.y * S[1];
    pPrime.z = p.z * S[2];
    return pPrime;
}

Point rotateOnX(Point p, float theta) {
    Point pPrime;
    pPrime.x = p.x;
    pPrime.y = p.y * cos(degToRad(theta)) - p.z * sin(degToRad(theta));
    pPrime.z = p.y * sin(degToRad(theta)) + p.z * cos(degToRad(theta));
    return pPrime;
}

Point rotateOnY(Point p, float theta) {
    Point pPrime;
    pPrime.x = p.x * cos(degToRad(theta)) + p.z * sin(degToRad(theta));
    pPrime.y = p.y;
    pPrime.z = -p.x * sin(degToRad(theta)) + p.z * cos(degToRad(theta));
    return pPrime;
}

Point rotateOnZ(Point p, float theta) {
    Point pPrime;
    float thetaInRad = degToRad(theta);
    pPrime.x = p.x * cos(thetaInRad) - p.y * sin(thetaInRad);
    pPrime.y = p.x * sin(thetaInRad) + p.y * cos(thetaInRad);
    pPrime.z = p.z;
    return pPrime;
}

float degToRad(float deg) {
    return deg * PI / 180.0;
}

void printPoint(Point p) {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << p.x << " " << p.y << " " << p.z << std::endl;
}
