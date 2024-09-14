
#include <stdio.h>
#include <stdlib.h>

#define MAX_SIZE 100

// Função para trocar duas linhas de uma matriz
void swapRows(double matriz[MAX_SIZE][MAX_SIZE], int row1, int row2,
              int tamanho) {
  for (int i = 0; i < tamanho; i++) {
    double temp = matriz[row1][i];
    matriz[row1][i] = matriz[row2][i];
    matriz[row2][i] = temp;
  }
}

// Função para calcular o determinante da matriz
double Det_Matriz(double matrix[MAX_SIZE][MAX_SIZE], int size) {
  double det = 1;
  int sign = 1;

  // Eliminação de Gauss para triangularizar a matriz
  for (int i = 0; i < size - 1; i++) {
    if (matrix[i][i] == 0) {
      int j = i + 1;
      while (j < size && matrix[j][i] == 0) {
        j++;
      }
      if (j == size) {
        return 0; // Se todas as entradas abaixo da diagonal principal são zero,
                  // o determinante é zero
      }
      swapRows(matrix, i, j, size);
      sign *= -1;
    }
    for (int j = i + 1; j < size; j++) {
      double factor = matrix[j][i] / matrix[i][i];
      for (int k = i; k < size; k++) {
        matrix[j][k] -= factor * matrix[i][k];
      }
    }
  }

  // Multiplicar os elementos da diagonal principal
  for (int i = 0; i < size; i++) {
    det *= matrix[i][i];
  }

  return det * sign;
}

// ok
void inicializarMatriz(int colunas, int linhas,
                       double matriz[colunas][linhas]) {

  for (int i = 0; i < linhas; i++) {
    for (int j = 0; j < colunas; j++) {
      if (matriz[i][j] != 0) {
        matriz[i][j] = 0;
      }
    }
  }
}



// ok
void PrintarMatriz(int colunas, int linhas, double matriz[colunas][linhas]) {

  printf("\n\n --- Imprimir a Matriz --- \n");

  printf("\n ");
  for (int i = 0; i < linhas; i++) {
    for (int j = 0; j < colunas; j++) {
      printf("%.0f ", matriz[i][j]);
    }
    printf("\n ");
  }
  printf("\n");
}

void PrintarMatrizLU(int tamanho, double matriz[tamanho][tamanho]) {

  printf("\n\n --- Imprimir a Matriz LU --- \n");

  printf("\n ");
  for (int i = 0; i < tamanho; i++) {
    for (int j = 0; j < tamanho; j++) {
      printf("%.2f ", matriz[i][j]);
    }
    printf("\n ");
  }
  printf("\n");
}

// ok
void PrintarMatrizCompleta(int colunas, int linhas,
                           double matriz[colunas][linhas],
                           double matrizTermsIndep[colunas][linhas]) {

  printf("\n --- Imprimir a Matriz --- \n");

  colunas++;
  int linhaTermosIndep = linhas;
  printf("\n ");
  for (int i = 0; i < linhas; i++) {
    for (int j = 0; j < colunas; j++) {
      if (j == colunas - 1) {
        printf("| %.0f", matrizTermsIndep[i][0]);
      } else {
        printf("%.0f ", matriz[i][j]);
      }
    }
    printf("\n ");
  }
  printf("\n");
}

// ok
void PreencherTermosIndep(int colunas, int linhas,
                          double matrizTermsIndep[colunas][linhas]) {

  printf("\n\n --- Preencher a Matriz Dos Termos Independentes --- \n\n");

  for (int i = 0; i < linhas; i++) {
    printf("Valor da posição [%d][%d] -> ", i, colunas - 1);
    scanf("%lf", &matrizTermsIndep[i][colunas - 1]);
  }

  PrintarMatriz(colunas, linhas, matrizTermsIndep);
}

// ok
void PreencherMatriz(int tamanho, double matriz[tamanho][tamanho]) {

  printf("\n\n --- Preencher a Matriz --- \n\n");

  for (int i = 0; i < tamanho; i++) {
    for (int j = 0; j < tamanho; j++) {
      printf("Valor da posição [%d][%d] -> ", i, j);
      scanf("%lf", &matriz[i][j]);
    }
  }
  PrintarMatriz(tamanho, tamanho, matriz);
}

void decomposicaoLU(int size, double L[MAX_SIZE][MAX_SIZE],
                    double U[MAX_SIZE][MAX_SIZE],
                    double A[MAX_SIZE][MAX_SIZE]) {

  // Inicializa a matriz L com zeros e a diagonal com uns
  for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        L[i][j] = 0.0; //tenta botar um else if aqui, com a condição de (i < j) pra colocar o numero como zero
        // for (int j = 0; j < size; j++) {
        //   if (i == j) {
        //     L[i][j] = 1;
        //   }else if (i < j) {
        //     L[i][j] = 0;
        //   }
        // }
        //assim ó
      }
  }

  // Inicializa a matriz U com zeros
  for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
          U[i][j] = 0;
      }
  }

  printf("\n\nMatriz L ->");
  PrintarMatrizLU(size, L);
  printf("\n\nMatriz U ->");
  PrintarMatrizLU(size, U);

            for (int j = 0; j < size; j++) {
                for (int i = j; i < size; i++) {
                    int max_row = i;
                    double max_val = U[i][j];
                    for (int k = i + 1; k < size; k++) {
                        if (U[k][j] > max_val) {
                            max_val = U[k][j];
                            max_row = k;
                        }
                    }
                    if (max_val == 0.0) {
                        printf("Erro: matriz singular\n");
                        exit(1);
                    }

                    swap_rows(U, i, max_row, size);
                    swap_rows(L, i, max_row, size);
                    if (i > 0) {
                        swap_rows(A, i, max_row, size);
                    }

                    U[i][j] = A[i][j];
                    for (int k = i + 1; k < size; k++) {
                        L[k][j] = A[k][j] / U[i][j];
                        for (int l = j + 1; l < size; l++) {
                            A[k][l] -= L[k][j] * U[i][l];
                        }
                    }
                }
            }

            for (int i = 0; i < size; i++) {
                L[i][i] = 1.0;
            }

  //-----------------------------
}

int main(void) {

  int tamanho;

  do {
    printf("\nTamanho da matriz -> ");
    scanf("%d", &tamanho);
    if (tamanho < 2) {
      printf("\nTamanho da matriz deve ser maior que 2");
    }
  } while (tamanho < 2 || tamanho > MAX_SIZE);

  // Matrizes --->
  double matriz[MAX_SIZE][MAX_SIZE];
  double matrizTermsIndep[1][MAX_SIZE];
  double L[MAX_SIZE][MAX_SIZE];
  double U[MAX_SIZE][MAX_SIZE];

  inicializarMatriz(tamanho, tamanho, matriz);
  inicializarMatriz(1, tamanho, matrizTermsIndep);

  PreencherMatriz(tamanho, matriz);
  PreencherTermosIndep(1, tamanho, matrizTermsIndep);

  PrintarMatrizCompleta(tamanho, tamanho, matriz, matrizTermsIndep);

  decomposicaoLU(tamanho, L, U, matriz);

  printf("\n\nMatriz L ->");
  PrintarMatrizLU(tamanho, L);
  printf("\n\nMatriz U ->");
  PrintarMatrizLU(tamanho, U);

  double valorDet = Det_Matriz(matriz, tamanho);

  printf("\nValor DetMatriz -> %.1lf ", valorDet);



  if (valorDet != 0) {
    printf("\nentrou");
    // calcular a LU
    decomposicaoLU(tamanho, L, U, matriz);
  } else {
    printf("\nA matriz tem multiplas soluções!!");
  }

  

  printf("\n\n  ---  Fim  ---  \n");
  return 0;
}

/*
 -----  LINKS ---->

 http://www.ene.unb.br/gaborges/recursos/programacao/gmatrix/gmatrixdoc.pdf

*/

// ok
//  Função para trocar duas linhas de uma matriz
// void swapRows(int tamanho, int row1, int row2, int matrix[tamanho][tamanho])
// {
//   for (int i = 0; i < tamanho; i++) {
//     double temp = matrix[row1][i];
//     matrix[row1][i] = matrix[row2][i];
//     matrix[row2][i] = temp;
//   }
// }

// // ok
// //  Função para calcular o determinante da matriz
// int DetMatriz(int tamanho, int matrix[tamanho][tamanho]) {

//   printf("\n --- Determinante da Matriz --- \n");

//   double det = 1;
//   int sign = 1;

//   // Eliminação de Gauss para triangularizar a matriz
//   for (int i = 0; i < tamanho -1 ; i++) {
//     if (matrix[i][i] == 0) {
//       int j = i + 1;
//       while (j < tamanho && matrix[j][i] == 0) {
//         j++;
//       }
//       if (j == tamanho) {
//         return 0; // Se todas as entradas abaixo da diagonal principal são
//         zero, o determinante é zero
//       }
//       swapRows(tamanho, i, j, matrix);
//       sign *= -1;
//     }
//     for (int j = i + 1; j < tamanho; j++) {
//       double factor = matrix[j][i] / matrix[i][i];
//       for (int k = i; k < tamanho; k++) {
//         matrix[j][k] -= factor * matrix[i][k];
//       }
//     }
//   }

//   // Multiplicar os elementos da diagonal principal
//   for (int i = 0; i < tamanho; i++) {
//     det *= matrix[i][i];
//   }

//   return det * sign;
// }