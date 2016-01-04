#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>
#include <GL/glew.h>
#include <GL/freeglut.h>

#define PI 3.14159265

struct Point {
    double x;
    double y;
    double z;
};

typedef std::pair<std::vector<std::string>, std::vector<double>> Transformation;
std::vector<Point> oPoints;
std::vector<Point> pPrimes;

void process();
arma::mat getCompositeMatrix(std::vector<Transformation> &transformations);
bool anyOriginalPointIsOrigin();
double degToRad(double deg);
void printPoint(Point p);

void drawScene();
void resize(int w, int h);
void setup();

std::vector<Point> transform(std::vector<Point> &points, arma::mat &T);
Point translate(Point p, double D[]);
Point scale(Point p, double S[]);
Point rotateOnX(Point p, double theta);
Point rotateOnY(Point p, double theta);
Point rotateOnZ(Point p, double theta);

int main(int argc, char* argv[]) {
    process();

    /* OpenGL related calls. */
    glutInit(&argc, argv);
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

        /* Read all the pairs of transformations. */
        for (int i = 0; i < t; i++) {
            Transformation t;
            std::string tName;
            bool translateToOrigin = false;
            std::cin >> tName;

            /* Check if a translation to the origin is needed. */
            if (tName == "s" || tName == "r") {
                if (!anyOriginalPointIsOrigin()) {
                    Transformation tOrigin;
                    tOrigin.first.push_back("t");
                    tOrigin.second.push_back(-oPoints[0].x);
                    tOrigin.second.push_back(-oPoints[0].y);
                    tOrigin.second.push_back(-oPoints[0].z);
                    translateToOrigin = true;
                    transformations.push_back(tOrigin);
                }
            }

            /* Start inserting the current transformation. */
            t.first.push_back(tName);

            /* Treat translation and scaling differently because of the
             * number of parameters */
            if (t.first[0] == "t" || t.first[0] == "s") {
                double c;
                std::cin >> c;
                t.second.push_back(c);
                std::cin >> c;
                t.second.push_back(c);
                std::cin >> c;
                t.second.push_back(c);
            } else if (t.first[0] == "r") {
                double c;
                std::string axis;
                std::cin >> axis;
                t.first.push_back(axis);
                std::cin >> c;
                t.second.push_back(c);
            }
            transformations.push_back(t);

            /* Check if a translation back from the origin is required. */
            if (translateToOrigin) {
                Transformation tBack;
                tBack.first.push_back("t");
                tBack.second.push_back(oPoints[0].x);
                tBack.second.push_back(oPoints[0].y);
                tBack.second.push_back(oPoints[0].z);
                transformations.push_back(tBack);
            }
        }

        /* Print the origin points. */
        for (auto point : oPoints) {
            printPoint(point);
        }
        T = getCompositeMatrix(transformations);

        /* Print the composite transformation matrix. */
        T.print();

        /* Apply the transformations (with the composite matrix)
         * to each point. */
        for (auto r : transform(oPoints, T)) {
            pPrimes.push_back(r);
        }

        /* Print the final points. */
        for (auto point : pPrimes) {
            printPoint(point);
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
    for (size_t i = transformations.size(); i-- > 0; ) {
        if (transformations[i].first[0] == "t") {
            B = { {1, 0, 0, transformations[i].second[0]},
                  {0, 1, 0, transformations[i].second[1]},
                  {0, 0, 1, transformations[i].second[2]},
                  {0, 0, 0, 1} };
        } else if (transformations[i].first[0] == "s") {
            B = { {transformations[i].second[0], 0, 0, 0},
                  {0, transformations[i].second[1], 0, 0},
                  {0, 0, transformations[i].second[2], 0},
                  {0, 0, 0, 1} };
        } else if (transformations[i].first[0] == "r") {
            double theta = degToRad(transformations[i].second[0]);
            double cosTheta = cos(theta), sinTheta = sin(theta);
            if (transformations[i].first[1] == "x") {
                B = { {1, 0, 0, 0},
                      {0, cosTheta, -sinTheta, 0},
                      {0, sinTheta, cosTheta, 0},
                      {0, 0, 0, 1} };
            } else if (transformations[i].first[1] == "y") {
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

/* Traverses de vector of original points to check if any is the origin. */
bool anyOriginalPointIsOrigin() {
    for (auto point : oPoints) {
        if (point.x == 0 && point.y == 0 && point.z == 0)
            return true;
    }
    return false;
}

void drawScene() {
    glClear(GL_COLOR_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(10.0, 10.0, 8.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

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
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
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

double degToRad(double deg) {
    return deg * PI / 180.0;
}

void printPoint(Point p) {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << p.x << " " << p.y << " " << p.z << std::endl;
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

Point translate(Point p, double D[]) {
    Point pPrime;
    pPrime.x = p.x + D[0];
    pPrime.y = p.y + D[1];
    pPrime.z = p.z + D[2];
    return pPrime;
}

Point scale(Point p, double S[]) {
    Point pPrime;
    pPrime.x = p.x * S[0];
    pPrime.y = p.y * S[1];
    pPrime.z = p.z * S[2];
    return pPrime;
}

Point rotateOnX(Point p, double theta) {
    Point pPrime;
    pPrime.x = p.x;
    pPrime.y = p.y * cos(degToRad(theta)) - p.z * sin(degToRad(theta));
    pPrime.z = p.y * sin(degToRad(theta)) + p.z * cos(degToRad(theta));
    return pPrime;
}

Point rotateOnY(Point p, double theta) {
    Point pPrime;
    pPrime.x = p.x * cos(degToRad(theta)) + p.z * sin(degToRad(theta));
    pPrime.y = p.y;
    pPrime.z = -p.x * sin(degToRad(theta)) + p.z * cos(degToRad(theta));
    return pPrime;
}

Point rotateOnZ(Point p, double theta) {
    Point pPrime;
    double thetaInRad = degToRad(theta);
    pPrime.x = p.x * cos(thetaInRad) - p.y * sin(thetaInRad);
    pPrime.y = p.x * sin(thetaInRad) + p.y * cos(thetaInRad);
    pPrime.z = p.z;
    return pPrime;
}
