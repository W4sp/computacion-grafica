/**
 * @filename:		lsystems3d.cpp
 * @last_updated:	03.01.2016 16:28:43
 * @author:			Cristóbal Leiva Aburto.
 * Basado en el libro de A. Lindenmayer "The Algorithmic Beauty of Plants"
 * Para compilar: g++ lsystems3d.cpp -o lsystems3d -lm
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <stack>

#define PI 				3.14159265

#define DEBUG			1
#define DIM				3
#define MAX_DESC 		256
#define DEFAULT_STEP	1
#define DEFAULT_ANGLE	45

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

/* Operaciones sobre matrices */

/* Multiplicación de matrices */
void mat_by_mat(double result[DIM][DIM], double A[DIM][DIM], double B[DIM][DIM])
{
	int i, j, k;
	double sum;
	for(i = 0; i < DIM; i++)
	{
		for(j = 0; j < DIM; j++)
		{
			sum = 0.0;
			for(k = 0; k < DIM; k++)
				sum += A[i][k]*B[k][j];
			result[i][j] = sum;
		}
	}
}

/* Guarda el resultado de una matriz en otra */
void assign_mat(double result[DIM][DIM], double M[DIM][DIM])
{
	int i, j;
	for(i = 0; i < DIM; i++)
		for(j = 0; j < DIM; j++)
			result[i][j] = M[i][j];
}

/* Guarda el resultado de un vector en otro */
void assign_vec(double result[DIM], double A[DIM])
{
	int i;
	for(i = 0; i < DIM; i++)
		result[i] = A[i];
}

/* Imprimir matriz */
void print_mat(double M[DIM][DIM])
{
	int i, j;
	printf("\n");
	for(i = 0; i < DIM; i++)
	{
		for(j = 0; j < DIM; j++)
			printf("%f ", M[i][j]);
		printf("\n");
	}
	printf("\n");
}

/* Imprimir vector */
void print_vec(double A[DIM])
{
	int i;
	printf("\n");
	for(i = 0; i < DIM; i++)
		printf("%f ", A[i]);
	printf("\n\n");
}

/* Multiplica matriz por vector */
void mat_by_vec(double result[DIM], double M[DIM][DIM], double V[DIM])
{
	int i, j;
	for(i = 0; i < DIM; i++)
	{
		result[i] = 0;
		for(j = 0; j < DIM; j++)
			result[i] += M[i][j]*V[j];
	}
}

/* Suma dos vectores */
void sum_vec(double result[DIM], double A[DIM], double B[DIM])
{
	int i;
	for(i = 0; i < DIM; i++)
		result[i] = A[i] + B[i];
}

/*
 * Matrices de transformación
 */

/* Rotar en torno a eje U (Z) */
void Ru_matrix(double R[DIM][DIM], double angle)
{
	double alfa = (angle*PI)/180.0;
	double mat[DIM][DIM] = {
		{cos(alfa), sin(alfa), 0},
		{sin(alfa*-1.0), cos(alfa), 0},
		{0, 0, 1}
	};
	assign_mat(R, mat);
}

/* Rotar en torno a eje L (Y) */
void Rl_matrix(double R[DIM][DIM], double angle)
{
	double alfa = (angle*PI)/180.0;
	double mat[DIM][DIM] = {
		{cos(alfa), 0, sin(alfa*-1.0)},
		{0, 1, 0},
		{sin(alfa), 0, cos(alfa)}
	};
	assign_mat(R, mat);
}

/* Rotar en torno a eje H (X) */
void Rh_matrix(double R[DIM][DIM], double angle)
{
	double alfa = (angle*PI)/180.0;
	double mat[DIM][DIM] = {
		{1, 0, 0},
		{0, cos(alfa), sin(alfa*-1.0)},
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
void get_argument(char *desc, int start, double *arg, int *jump)
{
	int i = start;
	char str_arg[MAX_DESC+1];
	
	if(desc[i+1] == '(')
	{
		i += 2;
		while(desc[i] != ')')
		{
			str_arg[i-start-2] = desc[i];
			i++;
		}
		
		/* atof: convierte string a float/double */
		*arg = atof(str_arg);
		*jump = i-start;
	}
	else
		/* Un salto de 0 significa que no hay argumento */
		*jump = 0;
}

void read_desc(char *desc, double *P)
{
	char action;
	double arg;
	int i, jump = 0, push = 0, pop = 0;
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
	
	/* Estado inicial */
	assign_mat(EstadoActual.T, T);
	assign_vec(EstadoActual.P, P0);

	for(i = 0; i < strlen(desc); i++)
	{
		/* Lee caracter y verifica si existe argumento */
		action = desc[i];
		get_argument(desc, i, &arg, &jump);
		
		/* Obtiene estado actual */
		assign_mat(T, EstadoActual.T);
		assign_vec(P0, EstadoActual.P);
		
		/* Los casos se describen en "L-systems: from the Theory to Visual Models of Plants"
		 * Apartado num. 5: The turtle interpretation of L-systems */
		switch(action)
		{
			case 'F':
				/* Si no hay argumento, entonces tomar valor por defecto */
				if(!jump) arg = DEFAULT_STEP;
				
				/* Tamaño del segmento a dibujar */
				L[0] = arg;
				mat_by_vec(P, T, L);
				printf("Dibujar segmento (%f, %f, %f), ", P0[0], P0[1], P0[2]);
				sum_vec(P0, P, P0);
				printf("(%f, %f, %f)\n", P0[0], P0[1], P0[2]);
				break;
			case '+':
				if(!jump) arg = DEFAULT_ANGLE;
				
				Ru_matrix(R, arg);
				/* Aplicar la transformación R a T: T*R */
				mat_by_mat(M, T, R);
				assign_mat(T, M);
				
				if(DEBUG) printf("Rotar hacia izquierda en torno a eje U.  Ru(%f)\n", arg);
				break;
			case '-':
				if(!jump) arg = DEFAULT_ANGLE;
				
				Ru_matrix(R, arg*-1.0);
				/* Aplicar la transformación R a T: T*R */
				mat_by_mat(M, T, R);
				assign_mat(T, M);
				
				if(DEBUG) printf("Rotar hacia derecha en torno a eje U. Ru(-%f)\n", arg);
				break;
			case '&':
				if(!jump) arg = DEFAULT_ANGLE;
				
				Rl_matrix(R, arg);
				/* Aplicar la transformación R a T: T*R */
				mat_by_mat(M, T, R);
				assign_mat(T, M);
				
				if(DEBUG) printf("Rotar hacia izquierda en torno a eje L. Rl(%f)\n", arg);
				break;
			case '^':
				if(!jump) arg = DEFAULT_ANGLE;
				
				Rl_matrix(R, arg*-1.0);
				/* Aplicar la transformación R a T: T*R */
				mat_by_mat(M, T, R);
				assign_mat(T, M);
				
				if(DEBUG) printf("Rotar hacia derecha en torno a eje L. Rl(-%f)\n", arg);
				break;
			case '\\':
				if(!jump) arg = DEFAULT_ANGLE;
				
				Rh_matrix(R, arg);
				/* Aplicar la transformación R a T: T*R */
				mat_by_mat(M, T, R);
				assign_mat(T, M);

				if(DEBUG) printf("Rotar hacia izquierda en torno a eje H. Rh(%f)\n", arg);
				break;
			case '/':
				if(!jump) arg = DEFAULT_ANGLE;
				
				Rh_matrix(R, arg*-1.0);
				/* Aplicar la transformación R a T: T*R */
				mat_by_mat(M, T, R);
				assign_mat(T, M);
				
				if(DEBUG) printf("Rotar hacia derecha en torno a eje H. Rh(-%f)\n", arg);
				break;
			case '[':
			
				/* Seteamos flag para hacer push del estado actual */
				push = 1;
			
				if(DEBUG) printf("Guardar el estado actual en la pila.\n");
				break;
			case ']':
			
				/* Seteamos flag para hacer pop de la pila y actualizar estado */
				pop = 1;
			
				if(DEBUG) printf("Obtener estado desde la pila y actualizarlo como estado actual.\n");
				break;
			default:
				break;
		}
		i += jump;
		
		/* Guardar estado actual */
		assign_mat(EstadoActual.T, T);
		assign_vec(EstadoActual.P, P0);

		/* Verificamos si debemos sacar el estado actual desde la pila */
		if(pop && !PilaEstados.empty())
		{
			EstadoActual = PilaEstados.top();			
			PilaEstados.pop();
			pop = 0;
		}
		/* Verificamos si debemos guardar estado actual en la pila */
		else if(push)
		{
			PilaEstados.push(EstadoActual);
			push = 0;
		}
	}
	
	/* Vaciar pila */
	while(!PilaEstados.empty())
	{
		PilaEstados.pop();
	}
}

int main()
{
	/* Descripción textual del L-system */
	char lsystem_desc[MAX_DESC+1] = {"F(2)[-F[-F]F]/(137.5)F(1.5)[-F]F"};
	
	/* Punto inicial */
	double P[DIM] = {0.0, 0.0, 0.0};
	
	read_desc(lsystem_desc, P);
	
	return 0;
}
