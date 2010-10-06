//vector operations, written by Danny bickson

#ifndef _VEC_OPS_H
#define _VEC_OPS_H

#include <math.h>

namespace gl{

    sdouble square_sum(const sdouble * x, int len){
        sdouble xs = 0.0;
        for(int i=0; i<len; i++) {
            xs += x[i] * x[i];
        }
        return  xs;
    }

   template<typename parray>      
   inline double rmse(const double * x1, const parray * x2, int len, double val, double & sum){
	    
             //assert(val>=0 && val <= 5);
             sum = 0;
             double ret = 0;
             for (int i=0; i< len; i++){
                ret = (x1[i] * x2[i]);
                //assert(ret>0);  //TODO
                sum+= ret;
             }
             
	    return pow(sum - val, 2);
        }

    inline double rmse(const double * x1, const double * x2, const double *x3, int len, double val, double &sum){
            if (x3 == NULL) //matrix	
                 return rmse(x1,x2,len,val,sum);

             //assert(val>=0 && val <= 5);
	     sum = 0;
             double ret = 0;
             for (int i=0; i< len; i++){
                ret = (x1[i] * x2[i] * x3[i]);
                //assert(ret>0);  //TODO
                sum+= ret;
             }
             
	    return pow(sum - val, 2);
        }


   inline double dot3(const double * x1, const double * x2, const double * x3, int len){
	     double sum = 0;
             double ret = 0;
             for (int i=0; i< len; i++){
                ret = (x1[i] * x2[i] * x3[i]);
                //assert(ret>0);  //TODO
                sum+= ret;
             }
            return sum;
        }

  template<typename parray>
  inline void dot2(const double * x1, const parray * x3, double * ret, int len){
             for (int i=0; i< len; i++){
                ret[i] = (x1[i] * x3[i]);
                //assert(ret>0);  //TODO
             }
        }


    sdouble norm2(sdouble* x, int len) {
        sdouble xs = square_sum(x,len);
        return  sqrt(xs);
    }
    sdouble norm_inf(const sdouble* x, int len) {
        sdouble max = -1e100;
        for (int i=0; i<len; i++){
		if (fabs(x[i]) > max)
			max = fabs(x[i]);
	}
        return max;
    }
    sdouble norm2(const sdouble* x, int len, sdouble *y, int len2) {
        assert(len==len2);
        sdouble norm =0;
        for (int i=0; i<len; i++)
           norm += ((x[i]-y[i])*(x[i]-y[i]));
        return  sqrt(norm);
    }
     
     sdouble norm1(const sdouble* x, int len) {
        sdouble xs = 0.0;
        for(int i=0; i<len; i++) {
            xs += fabs(x[i]);
        }
        return  xs;
    }
     sdouble norm0(const sdouble* x, int len) {
        sdouble xs = 0.0;
        for(int i=0; i<len; i++) {
            xs += (x[i] != 0);
        }
        return  xs;
    }
     
    sdouble *clone(const sdouble *x, int len){
       sdouble * ret = new sdouble[len];
       memcpy( ret, x, sizeof(sdouble)*len);
       return ret;
    }
     
     sdouble arraymul(const sdouble* x, int lenx, sdouble* y, int leny) {
       assert(lenx==leny);
        sdouble sum = 0;
        for(int i=0; i< lenx; i++) 
        	sum += (x[i] * y[i]);
        return sum;
    }
    sdouble * vpow(sdouble* x, int len, int powl) {
        for(int i=0; i<len; i++) 
        	x[i] = pow(x[i], powl);
	return x;
    } 
     void arraymul(sdouble* x, int len, sdouble y) {
        for(int i=0; i<len; i++) 
        	x[i] *= y;
    }
    
     sdouble sum(sdouble* x, int len) {
        sdouble sum = 0;
        for(int i=0; i< len; i++) 
        	sum += x[i];
        return sum;
    }
    
     sdouble sum_log(sdouble * x, int len){
        sdouble sum = 0;
        for(int i=0; i< len; i++) 
        	sum += log(x[i]);
        return sum;


     }

     sdouble* minus(sdouble* x, int len) {
        //sdouble *ret = new sdouble[len];
    	for(int i=0; i<len; i++) 
        	x[i] = -x[i];
        return x;
    }

     sdouble* vlog(sdouble* x, int len) {
        for(int i=0; i< len; i++) 
        	x[i] = log(x[i]);
        return x;
    }

     sdouble max(sdouble* x, int len) {
        sdouble max = -1e100;
    	for(int i=0; i< len; i++) 
        	if (x[i] > max)
        		max = x[i];
        return max;
    }

    
        sdouble* concat(sdouble* dx, int lenx, sdouble* du, int lenu) {
             sdouble * ret = new sdouble[lenx+lenu];
	     memcpy(ret, dx, lenx*sizeof(sdouble));
	     memcpy(ret + lenx*sizeof(sdouble), du, lenu*sizeof(sdouble));
	     return ret;
	}


	 sdouble* mul(sdouble* dx, int len, sdouble s) {
		assert(s!= 0);//, "multiplication by zero");
		sdouble * ret = new sdouble[len];
		for (int i=0; i< len; i++)
			ret[i] = dx[i] *s;
		return ret;
	}
	 double* dot(const double* dx1, const double * dx2, int len) {
		double * ret = new double[len];
		for (int i=0; i< len; i++)
			ret[i] = dx1[i] *dx2[i];
		return ret;
	}



	 sdouble* add(sdouble* x, int lenx, sdouble* y, int leny) {
		assert(lenx == leny); 
		sdouble * ret = new sdouble[lenx];
		for (int i=0; i< lenx; i++)
			ret[i] = x[i] + y[i];
		return ret;
	}
	 sdouble * add(sdouble* x, int lenx, sdouble y) {
    	        for (int i=0; i< lenx; i++)
			x[i] += y;
                return x;
         }
	 sdouble* sub(sdouble* x, int lenx, sdouble* y, int leny) {
		assert(lenx == leny); 
		sdouble * ret = new sdouble[lenx];
		for (int i=0; i< lenx; i++)
			ret[i] = x[i] - y[i];
		return ret;
	}
        sdouble* _sub(sdouble* x, int lenx, sdouble* y, int leny) {
		assert(lenx == leny); 
		for (int i=0; i< lenx; i++)
			x[i] -= y[i];
		return x;
	}

	 sdouble* arraymulalloc(sdouble* nu, int lennu, sdouble d) {
		sdouble * ret = new sdouble[lennu];
		for (int i=0; i< lennu; i++)
			ret[i] = nu[i]*d;
		return ret;
	}


	 sdouble * div(sdouble* x, int lenx, sdouble* y, int leny) {
		assert(lenx == leny); //, "wrong len");
		sdouble * ret = new sdouble[lenx];
		for (int i=0; i< lenx; i++){
			assert(y[i] != 0);// "division by zero");
			ret[i] = x[i]/y[i];
		}
		return ret;
	}


	 sdouble vmin(sdouble* x, int len) {
		sdouble min = 1e100;
		for (int i=0; i< len; i++){
			if (min > x[i])
				min = x[i];
		}
		return min;
	}


	 sdouble* vsqrt(sdouble* x, int len) {
		sdouble * ret = new sdouble[len];
		for (int i=0; i< len; i++){
			assert(x[i] >= 0); //, "negative square root");
			ret[i] = sqrt(x[i]);
		}
		return ret;
	}


	 sdouble* inv(sdouble* x, int len) {
		sdouble * ret = new sdouble[len];
		for (int i=0; i< len; i++)
			ret[i] = 1.0/x[i];
		return ret;
	}


	 sdouble* square(sdouble* x, int len) {
		sdouble * ret = new sdouble[len];
		for (int i=0; i< len; i++)
			ret[i] = x[i]*x[i];
		return ret;
	}


	 sdouble* ones(int n) {
		assert(n>0);
		sdouble * ret = new sdouble[n];
		for (int i=0; i< n; i++)
			ret[i] = 1.0;
		return ret;
	}

	 sdouble* ones(int n, sdouble val) {
		assert(n>0);
		sdouble * ret = new sdouble[n];
		for (int i=0; i< n; i++)
			ret[i] = val;
		return ret;
	}
        void ones(double * x, int len){
             for (int i=0; i<len; i++)
                 x[i] = 1.0;
        }
        void ones(double * x, int len, double val){
             for (int i=0; i<len; i++)
                 x[i] = val;
        }

        void rand(double * x, int len){
	     for (int i=0; i<len; i++)
                 x[i] = drand48();
        }
       double rmse(double * x1, double * x2, double * x3, int len, double val){

	     double sum = 0;
             for (int i=0; i< len; i++)
                sum += (x1[i] * x2[i] * x3[i]);

	     sum = pow(sum - val, 2);
            return sum;
        }

	 sdouble* zeros(int n) {
		assert(n>0);
                sdouble * ret = new sdouble[n];
		for (int i=0; i< n; i++)
			ret[i] = 0.0;
		return ret;
	}

      

	 void reset(sdouble* x, int start, int end) {
		for (int i=start; i< end; i++)
			x[i] = 0;
	}


        sdouble median(sdouble x[], int n){
            int i,j;
            sdouble temp;
            for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
		   if (x[j]<x[i]){
                        temp=x[i];
		        x[i]=x[j];
                        x[j]=temp;
                   }

                }
           }
           if (n%2==0)
		return (x[n/2]+x[n/2-1])/2.0;
           else return x[n/2];
      }



} //end namespace gl
#endif
