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

//deallocate space used by the matrix.
void mat_free(MAT *m);


//matrix displacement.

//return a matrix pointing to the specific row of given matrix.
MAT mat_row(MAT m, size_t row);

//copies the the in destination. size must be same.  
void mat_copy(MAT dst, MAT src);


//scale the matrix by given factor;
void mat_scale(MAT dst, double scalar);

//return equivalent upper traingular matrix of given matrix.
void UTM(MAT dst, MAT src, unsigned char print_steps);

//returns the determinant of a squre matrix... 'NOTE:ASSERTION ENABLED IF MATRIX IS NOT A SQURE MATRIX.
double mat_det(MAT matrix);
//return equivalent lower traingular matrix of given matrix.
//void LTM(MAT dst, MAT src);



//Matrix operations.
//sum of two matrices.
void mat_sum(MAT dst, MAT src);

//multiplication of two matrices.
void mat_dot(MAT dst, MAT src1, MAT src2);

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
void UTM(MAT dst, MAT src, unsigned char print_steps){
  MATH_ASSERT(src.rows == src.cols); // given matrix must be the squre matrix.
  
  MATH_ASSERT(src.rows == dst.rows);  //Dimensions of source must be equal to dimensions of Destination.
  MATH_ASSERT(src.cols == dst.cols);
  mat_copy(dst, src);

  MAT rx = mat_alloc(1, dst.cols);
  double f;
  for(size_t i = 0 ; i < dst.rows-1; i++){
    for(size_t j = i+1; j < dst.rows; j++){
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
double mat_det(MAT matrix){
  double det = 1.0f;
  MAT utm = mat_alloc(matrix.rows , matrix.cols);
  UTM(utm, matrix, 1);
  for(size_t i = 0 ; i < matrix.rows; i++){
    det *= MAT_AT(utm,i,i);
  }
  return det;
}
//void LTM(MAT dst, MAT src){}


void mat_print(MAT mat, char *name){
   MATH_PRINTF("%s[\n",name);
   for(size_t i = 0 ; i < mat.rows; i++){
     for(size_t j = 0 ; j < mat.cols; j++){
        MATH_PRINTF("\t  %lf",MAT_AT(mat, i, j));       
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


#endif //MATHFLOW_IMPL
