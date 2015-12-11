#include <iostream>
#include <cmath>

#define PI 3.14159265

struct Point {
    float x;
    float y;
    float z;
};

Point translate(Point p, float D[]);
Point scale(Point p, float S[]);
Point rotateOnX(Point p, float theta);
Point rotateOnY(Point p, float theta);
Point rotateOnZ(Point p, float theta);
Point rotateOnAny(Point p, float theta);
float degToRad(float deg);
void printPoint(Point p);

int main() {
    Point p;
    float D[3] = {1, 1, 1}, S[3] = {2, 2, 2};
    p.x = 1;
    p.y = 1;
    p.z = 1;
    printPoint(translate(p, D));
    printPoint(scale(p, S));
    printPoint(rotateOnX(p, 90));
    return EXIT_SUCCESS;
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
    pPrime.x = p.x * cos(degToRad(theta)) - p.y * sin(degToRad(theta));
    pPrime.y = p.x * sin(degToRad(theta)) + p.y * cos(degToRad(theta));
    pPrime.z = p.z;
    return pPrime;
}

float degToRad(float deg) {
    return deg * PI / 180.0;
}

void printPoint(Point p) {
    std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
}
