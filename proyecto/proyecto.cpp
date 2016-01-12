/**
 * L-systems
 * Basado en el libro de A. Lindenmayer "The Algorithmic Beauty of Plants"
 *
 * Para compilar: make
 * Para ejecutar: ./proyecto < data/[0-8].txt
 */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <stack>
#include <vector>
#include <array>
#include <GL/glew.h>
#include <GL/freeglut.h>

#define PI              3.14159265
#define ESC             27
#define DEBUG           1
#define DIM             3
#define MAX_DESC        16500
#define DEFAULT_STEP    1
#define DEFAULT_ANGLE   45

/* 
 * Estructura para guardar estado actual del L-system.
 * Esta estructura consiste en:
 * T: matriz de transformación actual
 * P: punto actual
 */
 
typedef struct {
    double T[DIM][DIM];
    double P[DIM];
} State;

State EstadoActual;

/* Pila para guardar estado actual al iniciar una nueva 'rama' (branch) */
std::stack<State> PilaEstados;
std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> lines;

float XAngle = 0.0;
float YAngle = 0.0;

/* Descripción textual del L-system */
char lsystem_desc[MAX_DESC + 1];
double langle = DEFAULT_ANGLE;
double lstep = DEFAULT_STEP;

/* Punto inicial */
double P[DIM] = {0.0, 0.0, 0.0};

void drawScene();
void resize(int w, int h);
void keyInput(unsigned char key, int x, int y);
void setup();
void mat_by_mat(double result[DIM][DIM], double A[DIM][DIM], double B[DIM][DIM]);
void assign_mat(double result[DIM][DIM], double M[DIM][DIM]);
void assign_vec(double result[DIM], double A[DIM]);
void print_mat(double M[DIM][DIM]);
void print_vec(double A[DIM]);
void mat_by_vec(double result[DIM], double M[DIM][DIM], double V[DIM]);
void sum_vec(double result[DIM], double A[DIM], double B[DIM]);
void Ru_matrix(double R[DIM][DIM], double angle);
void Rl_matrix(double R[DIM][DIM], double angle);
void Rh_matrix(double R[DIM][DIM], double angle);
void get_argument(char *desc, int start, double *arg, int *jump);
void read_desc(char *desc, double *P);

void drawScene() {
    float distance = 15.0;
    float XRad = XAngle / 180 * PI;
    float YRad = YAngle / 180 * PI;
    float x = sin(XRad) * distance;
    float y = sin(YRad) * distance;
    float z = cos(XRad) * distance;

    float lightPos0[] = {3.0, 17.0, 5.0, 1.0};
    float matAmbAndDif1[] = {0.9, 0.0, 0.0, 1.0};
    float matAmbAndDif2[] = {0.0, 0.9, 0.0, 1.0};
    float matSpec[] = {1.0, 1.0, 1.0, 1.0};
    float matShine[] = {50.0};

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(x + 5.0, 20.0, z + 3.0,
              0.0, 15.0, 0.0,
              0.0, 1.0, 0.0);

    glDisable(GL_LIGHTING);

    glPushMatrix();
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
    glTranslatef(lightPos0[0], lightPos0[1], lightPos0[2]);
    glColor3f(1.0, 1.0, 1.0);
    glutWireSphere(0.05, 8, 8);
    glPopMatrix();

    glEnable(GL_LIGHTING);

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbAndDif1);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matSpec);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, matShine);

    glLineWidth(1.0f);

    /* Draw red lines to depict the axes. */
    glColor4f(1.0, 0.0, 0.0, 0.35);
    glBegin(GL_LINES);
        glVertex4f(0.0, 0.0, 0.0, 1);
        glVertex4f(1000.0, 0.0, 0.0, 1);
        glVertex4f(0.0, 0.0, 0.0, 1);
        glVertex4f(0.0, 1000.0, 0.0, 1);
        glVertex4f(0.0, 0.0, 0.0, 1);
        glVertex4f(0.0, 0.0, 1000.0, 1);
    glEnd();

    /* Draw blue lines to depict the negative part of the axes. */
    glColor4f(0.0, 0.0, 1.0, 0.35);
    glBegin(GL_LINES);
        glVertex4f(0.0, 0.0, 0.0, 1);
        glVertex4f(-1000.0, 0.0, 0.0, 1);
        glVertex4f(0.0, 0.0, 0.0, 1);
        glVertex4f(0.0, -1000.0, 0.0, 1);
        glVertex4f(0.0, 0.0, 0.0, 1);
        glVertex4f(0.0, 0.0, -1000.0, 1);
    glEnd();

    glColor4f(0.0, 0.0, 1.0, 0.35);
    glBegin(GL_POLYGON);
        glVertex4f(0.0, 0.0, -10.0, 1);
        glVertex4f(10.0, 0.0, -10.0, 1);
        glVertex4f(10.0, 10.0, -10.0, 1);
        glVertex4f(0.0, 10.0, -10.0, 1);
    glEnd();

    /* Se renderizan los segmentos que conforman el fractal. */
    glColor4f(0.0, 1.0, 1.0, 1.0);
    glBegin(GL_LINES);
        for (auto l : lines) {
            glVertex4f(l.first[0], l.first[1], l.first[2], 1.0);
            glVertex4f(l.second[0], l.second[1], l.second[2], 1.0);
        }
    glEnd();

    glTranslatef(5.0, 15.0, 3.0);
    glutSolidSphere(1.5, 200, 200);

    glutSwapBuffers();
}

void resize(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-10.0, 10.0, -10.0, 10.0, 10.0, 50.0);
    glMatrixMode(GL_MODELVIEW);
}

void keyInput(unsigned char key, int x, int y) {
    switch (key) {
        case ESC:
            exit(EXIT_SUCCESS);
            break;
        default:
            break;
    }
}

void specialKeyInput(int key, int x, int y) {
    if (key == GLUT_KEY_UP) YAngle += 5;
    if (key == GLUT_KEY_DOWN) YAngle -= 5;
    if (key == GLUT_KEY_LEFT) XAngle -= 5;
    if (key == GLUT_KEY_RIGHT) XAngle += 5;
    glutPostRedisplay();
}

void setup() {
    /* Color de fondo negro. */
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_DEPTH_TEST);

    /* Activar iluminación de OpenGL. */
    glEnable(GL_LIGHTING);

    float lightAmb[] = {0.0, 0.0, 0.0, 1.0};
    /* Usualmente diffuse y specular se asignan a una luz brillante
     * para que no se altere el color de los objetos. */
    float lightDifAndSpec[] = {1.0, 1.0, 1.0, 1.0};
    float globAmb[] = {0.2, 0.2, 0.2, 1.0};

    /* Se setean las propiedades de luz de la fuente de luz 0. */
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDifAndSpec);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightDifAndSpec);

    glEnable(GL_LIGHT0); // Activar luz 0.
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globAmb);
}

int main(int argc, char* argv[]) {

    /* Leer parámetros del L-system desde stdin */
    scanf("%lf", &lstep);
    scanf("%lf", &langle);
    std::cin >> lsystem_desc;

    read_desc(lsystem_desc, P);

    /* OpenGL related calls. */
    glutInit(&argc, argv);
    /*glutInitContextVersion(2, 1);
    glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);*/
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Proyecto de Computación Gráfica");

    glutDisplayFunc(drawScene);
    glutReshapeFunc(resize);
    glutKeyboardFunc(keyInput);
    glutSpecialFunc(specialKeyInput);

    glewInit();

    setup();

    glutMainLoop();

    return EXIT_SUCCESS;
}

/* Operaciones sobre matrices */

/* Multiplicación de matrices */
void mat_by_mat(double result[DIM][DIM], double A[DIM][DIM], double B[DIM][DIM]) {
    int i, j, k;
    double sum;
    for (i = 0; i < DIM; i++) {
        for (j = 0; j < DIM; j++) {
            sum = 0.0;
            for (k = 0; k < DIM; k++) {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }
}

/* Guarda el resultado de una matriz en otra */
void assign_mat(double result[DIM][DIM], double M[DIM][DIM]) {
    int i, j;
    for (i = 0; i < DIM; i++)
        for (j = 0; j < DIM; j++)
            result[i][j] = M[i][j];
}

/* Guarda el resultado de un vector en otro */
void assign_vec(double result[DIM], double A[DIM]) {
    int i;
    for (i = 0; i < DIM; i++)
        result[i] = A[i];
}

/* Imprimir matriz */
void print_mat(double M[DIM][DIM]) {
    int i, j;
    printf("\n");
    for (i = 0; i < DIM; i++) {
        for (j = 0; j < DIM; j++)
            printf("%f ", M[i][j]);
        printf("\n");
    }
    printf("\n");
}

/* Imprimir vector */
void print_vec(double A[DIM]) {
    int i;
    printf("\n");
    for (i = 0; i < DIM; i++)
        printf("%f ", A[i]);
    printf("\n\n");
}

/* Multiplica matriz por vector */
void mat_by_vec(double result[DIM], double M[DIM][DIM], double V[DIM]) {
    int i, j;
    for (i = 0; i < DIM; i++) {
        result[i] = 0;
        for (j = 0; j < DIM; j++)
            result[i] += M[i][j] * V[j];
    }
}

/* Suma dos vectores */
void sum_vec(double result[DIM], double A[DIM], double B[DIM]) {
    int i;
    for (i = 0; i < DIM; i++)
        result[i] = A[i] + B[i];
}

/*
 * Matrices de transformación
 */

/* Rotar en torno a eje U (Z) */
void Ru_matrix(double R[DIM][DIM], double angle) {
    double alfa = (angle * PI)/180.0;
    double mat[DIM][DIM] = {
        {cos(alfa), sin(alfa), 0},
        {sin(alfa * -1.0), cos(alfa), 0},
        {0, 0, 1}
    };
    assign_mat(R, mat);
}

/* Rotar en torno a eje L (Y) */
void Rl_matrix(double R[DIM][DIM], double angle) {
    double alfa = (angle * PI)/180.0;
    double mat[DIM][DIM] = {
        {cos(alfa), 0, sin(alfa * -1.0)},
        {0, 1, 0},
        {sin(alfa), 0, cos(alfa)}
    };
    assign_mat(R, mat);
}

/* Rotar en torno a eje H (X) */
void Rh_matrix(double R[DIM][DIM], double angle) {
    double alfa = (angle * PI)/180.0;
    double mat[DIM][DIM] = {
        {1, 0, 0},
        {0, cos(alfa), sin(alfa * -1.0)},
        {0, sin(alfa), cos(alfa)}
    };
    assign_mat(R, mat);
}

/*
 * Leer el string de descripción desde el índice 'start' hasta encontrar
 * un cierre de paréntesis.  Asumimos que la string es una cadena bien
 * formada (con paréntesis balanceados).
 * Convertir el argumento leído entre paréntesis a un número flotante
 * (double). Además, calcular la cantidad de caractéres que hemos leído.
 * Si no hay argumento, 0 caracteres son leídos.
 */
void get_argument(char *desc, int start, double *arg, int *jump) {
    int i = start;
    char str_arg[MAX_DESC+1];

    if (desc[i+1] == '(') {
        i += 2;
        while (desc[i] != ')') {
            str_arg[i - start - 2] = desc[i];
            i++;
        }

        /* atof: convierte string a float/double */
        *arg = atof(str_arg);
        *jump = i - start;
    }
    else
        /* Un salto de 0 significa que no hay argumento */
        *jump = 0;
}

void read_desc(char *desc, double *P) {
    char action;
    double arg;
    int jump = 0, push = 0, pop = 0;
    /* Matriz para almacenar transformaciones a lo largo de iteraciones
     * En un comienzo apunta hacia Y+, ya que la primera columna indica
     * el Heading (hacia dónde apunta), la seguna cuál es la dirección 
     * hacia la izquierda (en este caso hacia X+) y cual es la dirección
     * hacia arriba (en este caso Z+). Cuando fue descrita por Lindenmayer
     * et al en "The Algorithmic Beauty of Plants" esta matriz es mencionada
     * como [H L U] (por Heading, Left, Up). */
    double T[DIM][DIM] = {{0, 1, 0}, {1, 0, 0}, {0, 0, 1}};
    /* Matriz para guardar rotaciones y resultados de operaciones */
    double R[DIM][DIM];
    double M[DIM][DIM];
    /* Vector para almacenar tamaño de segmento a dibujar */
    double L[DIM] = {0, 0, 0};
    /* Punto inicial */
    double P0[DIM] = {P[0], P[1], P[2]};

    std::pair<std::array<double, 3>, std::array<double, 3>> line;

    /* Estado inicial */
    assign_mat(EstadoActual.T, T);
    assign_vec(EstadoActual.P, P0);

    for (unsigned int i = 0; i < strlen(desc); i++) {
        /* Lee caracter y verifica si existe argumento */
        action = desc[i];
        get_argument(desc, i, &arg, &jump);

        /* Obtiene estado actual */
        assign_mat(T, EstadoActual.T);
        assign_vec(P0, EstadoActual.P);

        /* Los casos se describen en "L-systems: from the Theory to Visual Models of Plants"
         * Apartado num. 5: The turtle interpretation of L-systems */
        switch (action) {
            case 'F':
                /* Si no hay argumento, entonces tomar valor por defecto */
                if (!jump) arg = lstep;

                /* Tamaño del segmento a dibujar */
                L[0] = arg;
                mat_by_vec(P, T, L);
                printf("Dibujar segmento (%f, %f, %f), ", P0[0], P0[1], P0[2]);
                line.first = {P0[0], P0[1], P0[2]};
                sum_vec(P0, P, P0);
                printf("(%f, %f, %f)\n", P0[0], P0[1], P0[2]);
                line.second = {P0[0], P0[1], P0[2]};
                lines.push_back(line);
                break;
            case '+':
                if (!jump) arg = langle;

                Ru_matrix(R, arg);
                /* Aplicar la transformación R a T: T*R */
                mat_by_mat(M, T, R);
                assign_mat(T, M);

                if (DEBUG) printf("Rotar hacia izquierda en torno a eje U.  Ru(%f)\n", arg);
                break;
            case '-':
                if (!jump) arg = langle;

                Ru_matrix(R, arg*-1.0);
                /* Aplicar la transformación R a T: T*R */
                mat_by_mat(M, T, R);
                assign_mat(T, M);

                if (DEBUG) printf("Rotar hacia derecha en torno a eje U. Ru(-%f)\n", arg);
                break;
            case '&':
                if (!jump) arg = langle;

                Rl_matrix(R, arg);
                /* Aplicar la transformación R a T: T*R */
                mat_by_mat(M, T, R);
                assign_mat(T, M);

                if (DEBUG) printf("Rotar hacia izquierda en torno a eje L. Rl(%f)\n", arg);
                break;
            case '^':
                if (!jump) arg = langle;

                Rl_matrix(R, arg*-1.0);
                /* Aplicar la transformación R a T: T*R */
                mat_by_mat(M, T, R);
                assign_mat(T, M);

                if (DEBUG) printf("Rotar hacia derecha en torno a eje L. Rl(-%f)\n", arg);
                break;
            case '\\':
                if (!jump) arg = langle;

                Rh_matrix(R, arg);
                /* Aplicar la transformación R a T: T*R */
                mat_by_mat(M, T, R);
                assign_mat(T, M);

                if (DEBUG) printf("Rotar hacia izquierda en torno a eje H. Rh(%f)\n", arg);
                break;
            case '/':
                if (!jump) arg = langle;

                Rh_matrix(R, arg*-1.0);
                /* Aplicar la transformación R a T: T*R */
                mat_by_mat(M, T, R);
                assign_mat(T, M);

                if (DEBUG) printf("Rotar hacia derecha en torno a eje H. Rh(-%f)\n", arg);
                break;
            case '[':

                /* Seteamos flag para hacer push del estado actual */
                push = 1;

                if (DEBUG) printf("Guardar el estado actual en la pila.\n");
                break;
            case ']':

                /* Seteamos flag para hacer pop de la pila y actualizar estado */
                pop = 1;

                if (DEBUG) printf("Obtener estado desde la pila y actualizarlo como estado actual.\n");
                break;
            default:
                break;
        }
        i += jump;

        /* Guardar estado actual */
        assign_mat(EstadoActual.T, T);
        assign_vec(EstadoActual.P, P0);

        /* Verificamos si debemos sacar el estado actual desde la pila */
        if (pop && !PilaEstados.empty()) {
            EstadoActual = PilaEstados.top();
            PilaEstados.pop();
            pop = 0;
        }
        /* Verificamos si debemos guardar estado actual en la pila */
        else if (push) {
            PilaEstados.push(EstadoActual);
            push = 0;
        }
    }

    /* Vaciar pila */
    while (!PilaEstados.empty()) {
        PilaEstados.pop();
    }
}
