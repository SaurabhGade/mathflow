/*
 * mathflow.h
 * 
 * Author: Saurabh Gade <gadesaurabh3@gmail.com>
 * Date: 2024-08-13
 * Description: This header file defines the MathFlow library for linear algebra operations in C.
 * 
 * The MathFlow library provides essential methods for matrix operations, determinants, 
 * inversion, eigenvalues, eigenvectors, solvers, ...
 * 
 * License: MIT License
 */

#ifndef MATHFLOW
#define MATHFLOW

#ifndef MATH_PRINTF
#include<stdio.h>
#define MATH_PRINTF printf
#define MATH_SNPRINTF snprintf
#define MATH_SCANF scanf
#endif //MATH_PRINTF

#ifndef MATH_MALLOC
#include<stdlib.h>
#define MATH_MALLOC malloc
#define MATH_CALLOC calloc

#endif //MATH_MALLOC

#ifndef MATH_ASSERT
#include<assert.h>
#define MATH_ASSERT assert
#endif //MATH_ASSERT

typedef unsigned char bool;
#define true 0
#define false 1

typedef struct{
  size_t rows;
  size_t cols;
  size_t strid; 
  double *es;
} MAT;

//Access element at i,j index of given matrix.
#define MAT_AT(m, i, j) (m).es[(m).strid*(i)+(j)]
//instance token for MAT should be the name of the matrix.
#define MAT_PRINT(m) mat_print(m , #m);


//extra utilitis.
double rand_double(void);



//Matrix allocation.
//only method which can allocate memory.
MAT mat_alloc(size_t rows, size_t cols);

//initialize matrix with given values
void mat_init(MAT m , double val);

//initilize matrix with random flot values with given range
void mat_randf(MAT m, double upper_bound, double lower_bound);

//initilize matrix with random integer values with given range
void mat_rand(MAT m, double upper_bound, double lower_bound);

//initilize matrix by std input. 
void mat_scan(MAT m);


//Returns the identity matrix of given dimensions.
MAT mat_identity(size_t dim);


//deallocate space used by the matrix.
void mat_free(MAT *m);






//matrix displacement.

//return a matrix pointing to the specific row of given matrix.
MAT mat_row(MAT m, size_t row);

//copies the the in destination. size must be same.  
void mat_copy(MAT dst, MAT src);

//update row the ith row from the matrix by the given row .
void mat_update_row(MAT dst, MAT row, size_t i);

// shuffle matrix such as element at (i,j) should not zero if possible. (show rows) 
bool mat_shuffle_row(MAT m, size_t i, size_t j);

// swap two rows from given matrix 
 void mat_swap_row(MAT m, size_t i, size_t j);






//Matrix operations.
//sum of two matrices. 'WARNING: ASSERTION ENABLED IF, BOTH MATRICES ARE NOT OF SAME DIMENSIONS.

void mat_sum(MAT dst, MAT src);

//multiplication of two matrices. 'WARNING: ASSERTION ENABLED IF, COLUMNS OF SRC1 MATRIX IS NOT EQUEL TO ROWS OF SRC2 MATRIX. 
void mat_dot(MAT dst, MAT src1, MAT src2);

//scale the matrix by given factor;
void mat_scale(MAT dst, double scalar);

//return equivalent upper traingular matrix of given matrix.'WARNING: ASSERTION ENABLED IF, MATRIX IS NOT A SQURE MATRIX.
void mat_utm(MAT dst, MAT src, bool print_steps);

//return equivalent lower traingular matrix of given matrix.'WARNING: ASSERTION ENABLED IF, MATRIX IS NOT A SQURE MATRIX.
void mat_ltm(MAT dst, MAT src, bool print_steps);

//returns the determinant of a squre matrix... 'WARNING: ASSERTION ENABLED IF, MATRIX IS NOT A SQURE MATRIX.
double mat_det(MAT matrix, bool print_steps);



//convert given matrix in it's transpose.'WARNING: ASSERTION ENABLED IF, COLUMNS OF SRC1 MATRIX IS NOT EQUEL TO ROWS OF SRC2 MATRIX. 
void mat_transpose(MAT m);




// 'TODO: system of linear equitions.  
void mat_sys_linear_eq(MAT ans, MAT lhs, MAT rhs);




//print matrix.
void mat_print(MAT mat, char* name);


#endif //MATHFLOW



#ifndef MATHFLOW_IMPL

double rand_double(void){
  return (double)rand()/(double)RAND_MAX; 
}




MAT mat_alloc(size_t rows, size_t cols){
  MAT m;
  m.rows = rows;
  m.cols = cols;
  m.strid = cols;
  m.es = MATH_CALLOC(rows*cols,sizeof(*m.es));
  MATH_ASSERT(m.es != NULL);
  return m;
}
void mat_init(MAT m , double val){
  for(size_t i = 0 ; i < m.rows; i++){ for(size_t j = 0 ; j < m.cols; j++){
      MAT_AT(m, i, j) = val;
    }
  }
}

void mat_scan(MAT m){
  MATH_ASSERT(m.es != NULL);
  MATH_PRINTF("Enter the matrix of size %zux%zu\n",m.rows, m.cols); 
  for(size_t i = 0 ; i < m.rows ; i++){
    for(size_t j = 0 ; j < m.cols; j++){
     // MATH_PRINTF("(%zu,%zu) = ",i,j);
      MATH_SCANF("%lf",&MAT_AT(m, i, j));
    }
  }
}


MAT mat_identity(size_t dim){
  MAT id = mat_alloc(dim, dim);
  for(size_t i = 0 ; i < dim ; i++){
    MAT_AT(id, i, i) = 1; 
  }
  return id;
}


void mat_free(MAT *m){
  if(m->es == NULL)
    return;
  m->rows = 0;
  m->cols = 0;
  free(m->es);
}
void mat_randf(MAT m, double upper_bound, double lower_bound){
  for(size_t i = 0 ; i < m.rows; i++){
    for(size_t j = 0 ; j < m.cols; j++){
      MAT_AT(m, i, j) = rand_double() * (upper_bound-lower_bound)+lower_bound;
    }
  }
}

void mat_rand(MAT m, double upper_bound, double lower_bound){
  for(size_t i = 0 ; i < m.rows; i++){
    for(size_t j = 0 ; j < m.cols; j++){
      MAT_AT(m, i, j) =(int) (rand_double() * (upper_bound-lower_bound)+lower_bound); 
    }
  }
}

void mat_sum(MAT dst, MAT src){
  MATH_ASSERT(dst.rows == src.rows);
  MATH_ASSERT(dst.cols == src.cols);
  for(size_t i = 0 ; i < dst.rows; i++){
    for(size_t j = 0; j < dst.cols; j++){
      MAT_AT(dst, i, j) += MAT_AT(src, i, j);
    }
  }
}
void mat_dot(MAT dst, MAT src1, MAT src2){
  MATH_ASSERT(src1.cols == src2.rows);
  MATH_ASSERT(dst.rows == src1.rows);
  MATH_ASSERT(dst.cols == src2.cols);
  size_t n = src1.cols; 
  for(size_t i = 0 ; i < dst.rows; i++){
    for(size_t j = 0 ; j < dst.cols; j++){
      MAT_AT(dst, i, j) = 0;
      for(size_t k = 0 ; k < n ; k++){
        MAT_AT(dst, i, j) += MAT_AT(src1, i, k) * MAT_AT(src2,k, j); 
      }
    }
  }   
}

void mat_scale(MAT dst, double scalar){
  for(size_t i = 0 ; i < dst.rows; i++){
    for(size_t j = 0 ; j < dst.cols; j++){
      MAT_AT(dst, i, j) *= scalar;
    }
  }
}

//BUG: needed to fix. STILL BUT NOT FIXED
void mat_utm(MAT dst, MAT src, bool print_steps){
  MATH_ASSERT(src.rows == src.cols); // given matrix must be the squre matrix.
  
  MATH_ASSERT(src.rows == dst.rows);  //Dimensions of source must be equal to dimensions of Destination.
  MATH_ASSERT(src.cols == dst.cols);
  mat_copy(dst, src);

  MAT rx = mat_alloc(1, dst.cols);
  double f;
  for(size_t i = 0 ; i < dst.rows-1; i++){
    if(MAT_AT(dst, i, i) == 0){ 
      if(mat_shuffle_row(dst, i, i)) MAT_PRINT(dst);
    }
    for(size_t j = i+1; j < dst.rows; j++){
      
     // if(!MAT_AT(dst, i, j)) continue;
      mat_copy(rx, mat_row(dst, i));
      if(MAT_AT(rx, 0, i) == 0)
      f = 0;
      else
      f = MAT_AT(dst, j, i)/MAT_AT(rx, 0, i);

      if( (MAT_AT(rx , 0, i) < 0 && MAT_AT(dst, j, i) >= 0) || (MAT_AT(rx, 0, i) >= 0 && MAT_AT(dst, j, i) < 0) ) {
        if(f < 0) f *= -1;
      }else{
        if(f > 0) f *= -1;
      }
      mat_scale(rx, f);
      mat_sum(mat_row(dst,j), rx);
      if((int)print_steps)
        MAT_PRINT(dst);
    }
  }
  mat_free(&rx);
}

void mat_ltm(MAT dst, MAT src, bool print_steps){
  MATH_ASSERT(src.rows == src.cols); // given matrix must be the squre matrix.
  
  MATH_ASSERT(src.rows == dst.rows);  //Dimensions of source must be equal to dimensions of Destination.
  MATH_ASSERT(src.cols == dst.cols);
   
  mat_copy(dst, src);
  double f;
  MAT rx = mat_alloc(1, dst.cols);
  for(long i = dst.rows-1; i >= 0 ; i--){
    for(long j = i-1; j >= 0; j--){
      mat_copy(rx, mat_row(dst, i));
      f = (MAT_AT(rx, 0, i))? MAT_AT(dst,j,i)/MAT_AT(rx, 0, i): 0;
      if( (MAT_AT(rx , 0, i) < 0 && MAT_AT(dst, j, i) >= 0) || (MAT_AT(rx, 0, i) >= 0 && MAT_AT(dst, j, i) < 0) ) {
        if(f < 0) f *= -1;
      }else{
        if(f > 0) f *= -1;
      }
      mat_scale(rx, f);
      mat_sum(mat_row(dst, j), rx);
      if((int)print_steps)
        MAT_PRINT(dst);
    }
  }
  mat_free(&rx);
}


double mat_det(MAT matrix, bool print_steps){
  double det = 1.0f;
  MAT utm = mat_alloc(matrix.rows , matrix.cols);
  mat_utm(utm, matrix, print_steps);
  for(size_t i = 0 ; i < matrix.rows; i++){
    det *= MAT_AT(utm,i,i);
  }
  return det;
}

void mat_transpose(MAT m){
  MATH_ASSERT(m.rows == m.cols);

  for(size_t i = 0 ; i < m.rows; i++){
    for(size_t j = i+1 ; j < m.cols; j++){
        double d = MAT_AT(m, i, j);
        MAT_AT(m, i, j) = MAT_AT(m, j, i);
        MAT_AT(m, j, i) = d;
    }
  } 
}


void mat_print(MAT mat, char *name){
   MATH_PRINTF("%s[\n",name);
   for(size_t i = 0 ; i < mat.rows; i++){
     for(size_t j = 0 ; j < mat.cols; j++){
        double ld = MAT_AT(mat, i, j);
        if(ld >= -(1e-10) && ld <= 1e-10) ld = 0.f;
        //MATH_PRINTF("\t  %lf",MAT_AT(mat, i, j));       
        MATH_PRINTF("\t  %lf",ld);      
     }
    MATH_PRINTF("\n");
   }
  MATH_PRINTF("]\n");
}
MAT mat_row(MAT m, size_t row){
  MATH_ASSERT(row < m.rows);

  return (MAT){
    .rows = 1,
    .cols = m.cols,
    .strid = m.strid,
    .es = &MAT_AT(m, row, 0)
  };
}
void mat_copy(MAT dst, MAT src){
  MATH_ASSERT(dst.rows == src.rows);
  MATH_ASSERT(dst.cols == src.cols);
  for(size_t i = 0 ; i < dst.rows; i++){
    for(size_t j = 0 ; j < src.cols; j++){
      MAT_AT(dst, i, j) = MAT_AT(src, i, j); 
    }
  }
}

void mat_update_row(MAT dst, MAT row, size_t i){
  MATH_ASSERT(row.rows == 1);
  MATH_ASSERT(row.cols == dst.cols);
  MATH_ASSERT(i < dst.rows);
  for(size_t j = 0; j < row.cols; j++){
    MAT_AT(dst, i, j) = MAT_AT(row, 0, j);
  }
}

bool mat_shuffle_row(MAT m, size_t i, size_t j){
  MATH_ASSERT(i < m.rows);
  MATH_ASSERT(j < m.cols);

  if(MAT_AT(m, i, j)) return true;

  for(size_t k = i+1; k < m.rows; k++){
    if(MAT_AT(m, k, j) != 0){
      mat_swap_row(m, i, k);
      return true;
    }
  }
  return false;
}


void mat_swap_row(MAT m, size_t i, size_t j){
  MATH_ASSERT(i < m.rows);
  MATH_ASSERT(j < m.rows);
  if(i == j) return;
  MAT r = mat_alloc(1, m.cols);
  mat_copy(r, mat_row(m, i));
  mat_update_row(m, mat_row(m, j), i);
  mat_update_row(m, r, j);
  mat_free(&r);
}

#endif //MATHFLOW_IMPL
