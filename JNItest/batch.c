#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "batch.h"
#include "utils.h"

float ph2pr[128];
const char XDATA_pad[] = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC";
const char YDATA_pad[] = "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";

void fill_batch(t_batch * b, int x, int y, float initial,char *XDATA , char *YDATA,int * size,int * probs,int load,int pad,int batch) {
  t_inits * init = b->init;
  t_bbase * read = b->read;
  t_bbase * hapl = b->hapl;
  t_probs * prob = b->prob;
  
  int xp  = px(x,y);
  int xbp = pbp(xp);
 
  int yp  = py(y);
  int ybp = pbp(yp);
  int i;
  for(int i=0;i<128;i++)
  {
  	ph2pr[i]=powf(10.f, -((float)i) / 10.f);
 	
  } 
 
  init->batch_bytes = calc_batch_size(x,y,PES);
  init->x_size      = xp;
  init->x_padded    = xp;
  init->x_bppadded  = xbp;
  init->y_size      = yp;
  init->y_padded    = ybp;
  int pr_count		= 0;
  int init_prob_pointer=0;
  int init_hapl_pointer=0;
  int* size_temp =NULL;
  for(int f=0;f<batch;f++){
	init_prob_pointer+=size[f*2];
	init_hapl_pointer+=size[f*2+1];		
  }

  size_temp=(int*)(size+batch*2);			
  for (int k=0;k<PIPE_DEPTH;k++) {    
    if(k<PIPE_DEPTH-pad){	 
	init->initials[k] =(ldexpf(1.f,120))/(float)(k+1);
	int*   base_quals=(int*) (init_prob_pointer*4+probs+pr_count);
	int*   ins_quals =(int*) (init_prob_pointer*4+probs+size_temp[k*2]+pr_count);
	int*   del_quals =(int*) (init_prob_pointer*4+probs+size_temp[k*2]*2+pr_count);
	int*   gcp_quals =(int*) (init_prob_pointer*4+probs+size_temp[k*2]*3+pr_count);
	pr_count+=size_temp[k*2]*4;

	
   for (int i=0;i<xbp;i++) {	 
      if (i < size_temp[k*2] ) {
        read[i].base[k] = XDATA[i+init_prob_pointer];
		
      }	
      else // padding:
      {
        read[i].base[k] = 'S';
      }
    }
   for (int i=0;i<ybp;i++) {
      if (i < size_temp[k*2+1])
        hapl[i].base[k] =  YDATA[i+init_prob_pointer];
      else
        hapl[i].base[k] = 'S';
    }
  for (int i=0;i<xp;i++) {
		int _i = ins_quals[i] & 127;  
		int _d = del_quals[i] & 127;
		int _q = base_quals[i] & 127;	    
		float  distm = ph2pr[_q];

      prob[(i)*PIPE_DEPTH+k].p[0].b = FloatToUint(0.1) ; //zeta
      prob[(i)*PIPE_DEPTH+k].p[1].b = FloatToUint(ph2pr[_d]) ; //eta
      prob[(i)*PIPE_DEPTH+k].p[2].b = FloatToUint(0.1); //upsilon
      prob[(i)*PIPE_DEPTH+k].p[3].b = FloatToUint(ph2pr[_i]); //delta
      prob[(i)*PIPE_DEPTH+k].p[4].b = FloatToUint(0.9); //beta
      prob[(i)*PIPE_DEPTH+k].p[5].b = FloatToUint(1.0f-ph2pr[( _i+_d)&127]); //alpha
      prob[(i)*PIPE_DEPTH+k].p[6].b = FloatToUint(distm/3);  //distm_diff
      prob[(i)*PIPE_DEPTH+k].p[7].b = FloatToUint(1- distm); //distm_simi
	
    }
	}
  } 
  return;
}

uint32_t FloatToUint(float n)
{
   return (uint32_t)(*(uint32_t*)&n);
}

float UintToFloat(uint32_t n)
{
   return (float)(*(float*)&n);
}
void calculate_mids(t_batch * batch, int pair, int r, int c, float * M, float * I, float * D) {
  
  int w = c+1;
  t_inits * init = batch->init;
  t_bbase * read = batch->read;
  t_bbase * hapl = batch->hapl;
  t_probs * prob = batch->prob;  
    
  // Set to zero and intial value in the X direction
  for (int j = 0; j < c+1; j++)
  {
    M[j] = 0.0;
    I[j] = 0.0;
    D[j] = init->initials[pair];
  }

  // Set to zero in Y direction
  for (int i = 1; i < r+1; i++)
  {
    M[i*w] = 0.0;
    I[i*w] = 0.0;
    D[i*w] = 0.0;
  }
  
  for (int i = 1; i < r+1; i++) {
    for (int j = 1; j < c+1; j++) {
      
      float distm_simi = prob[(i-1)*PIPE_DEPTH+pair].p[7].f;
      float distm_diff = prob[(i-1)*PIPE_DEPTH+pair].p[6].f;
      float alpha      = prob[(i-1)*PIPE_DEPTH+pair].p[5].f;
      float beta       = prob[(i-1)*PIPE_DEPTH+pair].p[4].f;
      float delta      = prob[(i-1)*PIPE_DEPTH+pair].p[3].f;
      float upsilon    = prob[(i-1)*PIPE_DEPTH+pair].p[2].f;
      float eta        = prob[(i-1)*PIPE_DEPTH+pair].p[1].f;
      float zeta       = prob[(i-1)*PIPE_DEPTH+pair].p[0].f;
       
      float distm;
      unsigned char rb = read[i-1].base[pair];
      unsigned char hb = hapl[j-1].base[pair];
      
      if (rb == hb || rb == 'N' || hb == 'N') {
        distm = distm_simi;
      }
      else {
        distm = distm_diff;
      }
      
      //printf("read:%d,hapl:%d,read:%c,hapl:%c,distm:%8X,",i,j,r,h,*(uint32_t*)&distm);

      M[i*w+j] = distm * (alpha * M[(i-1)*w+(j-1)] + beta    * I[(i-1)*w+(j-1)] + beta * D[(i-1)*w+(j-1)]);
      I[i*w+j] =         (delta * M[(i-1)*w+(j  )] + upsilon * I[(i-1)*w+(j  )]                          );
      D[i*w+j] =         (eta   * M[(i  )*w+(j-1)]                              + zeta * D[(i  )*w+(j-1)]);
      
      //printf("m:%8X,i:%8X,d:%8X\n",*(uint32_t*)&M[i*w+j],*(uint32_t*)&I[i*w+j],*(uint32_t*)&D[i*w+j]);
      
      //printf("i:%2d,j:%2d,%8X %8X %8X\n",i,j,*(uint32_t*)&M[i*w+j],*(uint32_t*)&I[i*w+j],*(uint32_t*)&D[i*w+j]);
    }
    //printf("\n");
  }

  return;
}

int count_errors (uint32_t * hr, uint32_t * sr, int num_batches) {
	//int value_error = 0;
	int total_errors = 0;
  uint32_t hw;
  uint32_t sw;
  float hwf;
  float swf;

  for (int i = 0; i < num_batches; i++)  {
    for (int j = 0; j < PIPE_DEPTH; j++) {
      sw = sr[i*4*PIPE_DEPTH+j*4]; 
      hw = hr[i*4*PIPE_DEPTH+j*4];
      swf = *(float*)&sr[i*4*PIPE_DEPTH+j*4];
      hwf = *(float*)&hr[i*4*PIPE_DEPTH+j*4];
      float err = fabs(swf/hwf);
      
      if ((err < ERR_LOWER) || (err > ERR_UPPER))
      {
          total_errors++;
          printf("Error in batch %d pair %d:\nSW: %8X\nHW: %8X DIFF:%16f\n",i,j,sw,hw,err);
          exit(-1);
          DEBUG_PRINT("B%3d, P%3d, SW: %8X, HW: %8X, SW:%f, HW:%f\n",
                      i,
                      j,
                      sr[i*4*PIPE_DEPTH+j*4],
                      hr[i*4*PIPE_DEPTH+j*4],
                      *(float*)&sr[i*4*PIPE_DEPTH+j*4],
                      *(float*)&hr[i*4*PIPE_DEPTH+j*4]
          );
      }
    }    
  }

  return total_errors;
}




void init_batch_address(t_batch * b, void * batch, int x, int y)
{
    b->init = (t_inits*)((uint64_t) batch);
    b->read = (t_bbase*)((uint64_t) b->init + CACHELINE_BYTES);
    b->hapl = (t_bbase*)((uint64_t) b->read + (uint64_t) pbp(px(x,y)) * (uint64_t)sizeof(t_bbase));
    b->prob = (t_probs*)((uint64_t) b->hapl + (uint64_t) pbp(py(y)) * (uint64_t)sizeof(t_bbase));
}

size_t calc_batch_size(int x, int y, int e)
{
  
  // Pad X and Y if necessary
  int xpadded = px(x,y);
  int ypadded = py(y);
  
  // Pad X and Y basepairs
  int ybps    = pbp(ypadded);
  int xbps    = pbp(xpadded);
           
  //DEBUG_PRINT("Pair size                : Y=%3d, X=%3d \n", (int)y, (int)x);
  //DEBUG_PRINT("Padded size              : Y=%3d, X=%3d \n", (int)ypadded, (int)xpadded);
  
  size_t init_size = CACHELINE_BYTES;
  
  size_t xbases = xbps * sizeof(t_bbase);
  size_t ybases = ybps * sizeof(t_bbase);
    
  //DEBUG_PRINT("Size of intial values    : %6d / %d\n", (int) init_size, (int) init_size / CACHELINE_BYTES);
  //DEBUG_PRINT("Size of haplotype string : %6d / %d\n", (int) ybases, (int) ybases / CACHELINE_BYTES);
  //DEBUG_PRINT("Size of read string      : %6d / %d\n", (int) xbases, (int) xbases / CACHELINE_BYTES); 
  
  size_t prob_size = (xpadded) * PIPE_DEPTH * sizeof(t_probs);
  
  DEBUG_PRINT("Size of probabilities    : %6d / %d\n", (int) prob_size, (int) prob_size / CACHELINE_BYTES);

  // Determine size of batch in bytes.
  size_t batch_size = init_size + xbases + ybases + prob_size;
  
  //DEBUG_PRINT("Size of batch            : %6d / %d\n", (int) batch_size, (int) batch_size / CACHELINE_BYTES);
  
  return batch_size;
}

