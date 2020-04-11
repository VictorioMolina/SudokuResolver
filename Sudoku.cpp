/* --------------------------------- SUDOKU --------------------------------- */
#include <ga/GASimpleGA.h>      // Algoritmo Genetico Simple
#include <ga/GA1DArrayGenome.h> // Genoma --> array de enteros (dim. 1) alelos
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

using namespace std;

struct plantilla
{
    int tam;
    int *fijo;
};

// Prototipos de funcion
void leerSudoku(struct plantilla *, char *);
void InicioSudoku(GAGenome &);
int CruceSudoku(const GAGenome &, const GAGenome &, GAGenome *, GAGenome *);
bool checkColumna(int[], int *, int);
int MutacionSudoku(GAGenome &, float);
float Objective(GAGenome &);
GABoolean Termina(GAGeneticAlgorithm &);


int main(int argc, char **argv)
{
    // Leemos el Sudoku del fichero
    char* nombreF = argv[1];

    struct plantilla S;
    leerSudoku(&S, nombreF);

    cout << nombreF << " - " << S.tam << "x" << S.tam << endl; // Titulo del programa

    // Declaramos variables para los parametros del GA y las inicializamos
    int ngen = 12000;
    int popsize = atoi(argv[2]);
    char* selectorType = argv[3];
    float pcross = atof(argv[4]);
    float pmut = atof(argv[5]);

    cout << "Parametros:    - Tamano poblacion: " << popsize << endl;
    cout << "               - Tipo de seleccion: " << selectorType << endl;
    cout << "               - Probabilidad cruce: " << pcross << endl;
    cout << "               - Probabilidad mutacion: " << pmut << endl << endl;

    // Conjunto enumerado de alelos --> valores posibles de cada gen del genoma
    GAAlleleSet<int> alelos;
    for (int i = 0; i < S.tam; i++)
        alelos.add(i);

    // Creamos el genoma y definimos operadores de inicio, cruce y mutacion
    GA1DArrayAlleleGenome<int> genome(S.tam * S.tam, alelos, Objective, &S);
    genome.initializer(InicioSudoku);
    genome.crossover(CruceSudoku);
    genome.mutator(MutacionSudoku);

    // Creamos el algoritmo genetico
    GASimpleGA ga(genome);

    // Inicializamos - minimizar funcion objetivo, tamaño poblacion, num generaciones,
    // pr. cruce y pr. mutacion, seleccion y le indicamos que evolucione.

    ga.minimaxi(-1);
    ga.populationSize(popsize);
    ga.nGenerations(ngen);
    ga.pCrossover(pcross);
    ga.pMutation(pmut);

    if(strcmp(selectorType,"GARouletteWheelSelector") == 0){
        GARouletteWheelSelector selector;
        ga.selector(selector);
    }
    else if (strcmp(selectorType,"GATournamentSelector") == 0){
         GATournamentSelector selector;
         ga.selector(selector);
    }

    ga.terminator(Termina);
    ga.evolve(1);

    // Imprimimos el mejor individuo que encuentra el GA y su valor fitness
    cout << "El GA encuentra la solucion ( " << ga.statistics().bestIndividual() << ") con valor fitness " << ga.statistics().minEver() << endl;
    cout << "Numero de generaciones: " << ga.statistics().generation() << endl;
    return 0;
}

// Funcion de lectura del Sudoku
void leerSudoku(struct plantilla *S, char *nombreF)
{
    ifstream f(nombreF);
    f >> S->tam;
    S->fijo = new int[S->tam * S->tam];
    for (int i = 0; i < S->tam * S->tam; i++)
        f >> S->fijo[i];
    f.close();
}

// Funcion de inicio de Sudoku
void InicioSudoku(GAGenome &g)
{
    GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &)g;

    struct plantilla *datos;
    datos = (struct plantilla *) genome.userData();

    int aux[datos->tam];

    for (int f = 0; f < datos->tam; f++)
    {
        for (int j = 0; j < datos->tam; j++)
            aux[j] = 0;

        for (int j = 1; j <= datos->tam; j++)
        {
            int v = GARandomInt(0, datos->tam - 1);
            while (aux[v] != 0)
                v = (v + 1) % datos->tam;
            aux[v] = j;
        }


        int i = 0;
        while (i < datos->tam)
        {
            while ((datos->fijo[(f * datos->tam) + i] == 0) && (i < datos->tam))
                i++;

            if (i < datos->tam)
            {
                bool encontrado = false;
                for (int j = 0; (j < datos->tam) && (!encontrado); j++)
                    if (aux[j] == datos->fijo[(f * datos->tam) + i])
                    {
                        encontrado = true;
                        aux[j] = aux[i];
                    }
                aux[i] = datos->fijo[(f * datos->tam) + i];
            }
            i++;
        }

        for (int c = 0; c < datos->tam; c++)
            genome.gene((f * datos->tam) + c, aux[c]);
    }
}

// Funcion de Cruce de Sudoku
int CruceSudoku(const GAGenome &p1, const GAGenome &p2, GAGenome *c1, GAGenome *c2)
{
    const GA1DArrayAlleleGenome<int> &m = (GA1DArrayAlleleGenome<int> &)p1;
    const GA1DArrayAlleleGenome<int> &p = (GA1DArrayAlleleGenome<int> &)p2;

    struct plantilla *datos = (struct plantilla *)m.userData();
    int n = 0;

    int punto1 = GARandomInt(0, m.length());
    while ((punto1 % datos->tam) != 0)
        punto1++;
    int punto2 = m.length() - punto1;

    if (c1)
    {
        GA1DArrayGenome<int> &h1 = (GA1DArrayAlleleGenome<int>&) *c1;
        h1.copy(m, 0, 0, punto1); // el metodo copy esta definido en la clase GA1DArrayGenome
        h1.copy(p, punto1, punto1, punto2);
        n++;
    }

    if (c2)
    {
        GA1DArrayGenome<int> &h2 = (GA1DArrayAlleleGenome<int>&) *c2;
        h2.copy(p, 0, 0, punto1);
        h2.copy(m, punto1, punto1, punto2);
        n++;
    }

    return n;
}

bool checkColumna(int col[], int *check, int tam)
{
    bool repe = false;

    for (int i = 0; i < tam; i++)
        check[i] = 0;

    for (int i = 0; i < tam; i++)
        check[col[i] - 1]++;
    for (int i = 0; i < tam; i++)
        if (check[i] > 1)
            repe = true;

    return repe;
}

// Funcion de Mutacion
int MutacionSudoku(GAGenome &g, float pmut)
{
    GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &)g;

    struct plantilla *datos;
    datos = (struct plantilla *)genome.userData();
    int nmut = 0;
    int aux;
    int fil;
    bool fila;

    int caux[datos->tam];
    int *checkC = new int[datos->tam];

    if (pmut <= 0.0)
        return 0;

    for (int f = 0; f < datos->tam; f++)
        for (int c = 0; c < datos->tam; c++)
            if (datos->fijo[(f * datos->tam) + c] == 0)
            {
                if (GAFlipCoin(pmut))
                {
                    if (GAFlipCoin(0.5))
                        fila = true;
                    else
                        fila = false;

                    if (!fila)
                    {

                        for (int j = 0; j < datos->tam; j++)
                            caux[j] = genome.gene((j * datos->tam) + c);
                        if (checkColumna(caux, checkC, datos->tam))
                        {
                            int v1 = GARandomInt(0, datos->tam - 1);
                            while (checkC[v1] <= 1)
                                v1 = (v1 + 1) % datos->tam;
                            v1++;
                            int v2 = GARandomInt(0, datos->tam - 1);
                            while (checkC[v2] != 0)
                                v2 = (v2 + 1) % datos->tam;
                            v2++;

                            bool encontrado = false;
                            for (int j = 0; j < datos->tam && !encontrado; j++)
                                if ((datos->fijo[j * (datos->tam) + c] == 0) && (genome.gene(j * (datos->tam) + c) == v1))
                                {
                                    encontrado = true;
                                    genome.gene((j * datos->tam) + c, v2);
                                    fil = j;
                                }

                            int col = (c + 1) % datos->tam;
                            while (genome.gene((fil * datos->tam) + col) != v2)
                                col = (col + 1) % datos->tam;
                            if (datos->fijo[(fil * datos->tam) + col] == 0)
                            {
                                nmut++;
                                genome.gene((fil * datos->tam) + col, v1);
                            }
                            else
                            {
                                genome.gene((fil * datos->tam) + c, v1);
                            }
                        }
                    }
                    else
                    {
                        int v1 = (c + 1) % datos->tam;
                        while ((datos->fijo[(f * datos->tam) + v1] != 0))
                            v1 = (v1 + 1) % datos->tam;
                        aux = genome.gene((f * datos->tam) + c);
                        genome.gene((f * datos->tam) + c, genome.gene((f * datos->tam) + v1));
                        genome.gene((f * datos->tam) + v1, aux);
                        nmut++;
                    }
                }
            }

    return nmut;
}

// Funcion objetivo
float Objective(GAGenome &g)
{
    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;
    float nrepeticiones = 0;

    struct plantilla *datos;
    datos = (struct plantilla *)genome.userData();

    // Numero de filas y de columnas
    int nfilas = datos->tam;
    int ncolumnas = nfilas; // El sudoku es cuadrado

    int aux[nfilas][ncolumnas] = {0};

    /*
        Represento el genoma unidimensional
        como una matriz de dos dimensiones
        para facilitar la implementación del
        algoritmo.
    */
    for(int i = 0; i < nfilas; i++){
        for(int j = 0; j < ncolumnas; j++){
            aux[i][j] = genome.gene(j+i*ncolumnas);
        }
    }

    // No hay repeticiones por filas

    // Repeticiones en la misma columnas
   vector<int> encontrados(datos->tam); // Vector de encontrados
   fill(encontrados.begin(), encontrados.end(), 0); // Inicializamos el vector a 0
   for(int j = 0; j < ncolumnas; j++){
        for(int i = 0; i < nfilas; i++){
            encontrados[aux[i][j] - 1]++;
             // En caso de que el número haya aparecido más de una vez...
            if(encontrados[aux[i][j] - 1] > 1){
                // Aumentamos el número de repeticiones
                nrepeticiones++;
            }
        }
        fill(encontrados.begin(), encontrados.end(), 0); // Reiniciamos el vector
    }

    // Repeticiones en la misma cuadricula
    int tcuadricula = datos->tam/sqrt(datos->tam);
    for(int i = 0; i < sqrt(datos->tam); i++){ // Para cada fila de cuadriculas
        for(int j = 0; j < sqrt(datos->tam); j++){ // Para cada columna de cuadriculas
            for(int fila = 0; fila < tcuadricula; fila++){ // Para cada fila dentro de la cuadricula
                for(int columna = 0; columna < tcuadricula; columna++){ // Para cada columna dentro de la cuadricula
                    encontrados[aux[fila+tcuadricula*i][columna+tcuadricula*j] - 1]++;
                    // En caso de que el número haya aparecido más de una vez...
                    if(encontrados[aux[fila+tcuadricula*i][columna+tcuadricula*j] - 1] > 1){
                        // Aumentamos el número de repeticiones
                        nrepeticiones++;
                    }
                }
            }
            fill(encontrados.begin(), encontrados.end(), 0); // Reiniciamos el vector
        }
    }


    return nrepeticiones;
}

// Funcion de terminacion
GABoolean Termina(GAGeneticAlgorithm &ga)
{
    /* El algoritmo termina cuando se haya alcanzado el número de
    generaciones límite o cuando hayamos encontrado
    una solución (fitness=0) */
    if ((ga.statistics().minEver()==0) ||
        (ga.statistics().generation()==ga.nGenerations())) return gaTrue;

    return gaFalse;
}
