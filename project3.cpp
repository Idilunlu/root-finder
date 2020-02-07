#include <stdio.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

class all{
	public:
		float funct(float*command,int arg,float x);       // function to construct polynomials
		void bisection(float*array,float hval,float lval,float tol,float argc);   // algorithm for bisection
		void secant(float*array,float hval,float lval,float tol,float argc);   // algorithm for secant
		void hybrid(float*array,float hval,float lval,float tol,float argc);   // hybrid algorithm
};

 
float all::funct(float*command,int arg,float x){ 
    float result=0;
    int n=arg-4;          // gets the number of coefficient inputs in command line arguments
    for(int i=0;i<n;i++){
    float mult;  
	mult=pow(x,n-i-1);       // getting the decreasing powers of x with loop for polynomial equation
    result+=command[i]*mult;    // multiplying appropriate power of x with the corresponding coefficient of polynomial
} 
	return result;    
}

void all::bisection(float*array,float hval,float lval,float tol,float argc){
	
	    int bisectioncounter=0;         // will be used for limiting the loop
	    float f_lval,f_hval,f_root,root;   // polynomial results of variables
	    
	    f_hval=funct(array,argc,hval);     // getting f values of inital guesses and root;
	    f_lval=funct(array,argc,lval);
	    
	    if(f_lval==0 && f_hval==0){         
				cout<<"Initial guesses are roots themselves. "<<endl<<endl;	//	break;
			}
			
		else if(f_lval==0 || f_hval==0 ){  
				cout<<"One of the initial guesses is root itself. "<<endl<<endl;
	
		        root=lval+(hval-lval)/2;       // calculating the middle point
		
		    if(f_hval==0){         // eliminating initial guess corresponding to root
			hval=root;             // needed for bisection method to properly work
	     	}
		    else if(f_lval==0){      
			lval=root;
		    }
		
		    bisectioncounter++;    // actually one iteration is already made outside of the while loop

	     	while(hval-lval>tol){      // copmaring the distance to exit loop
		   
		    root=lval+(hval-lval)/2;  
		    
	         bisectioncounter++;
	    
	         	if(signbit(funct(array,argc,lval))==signbit(funct(array,argc,root))){     
			    lval=root;           // sign of value of border variables cant be same 
	        	}                    // new borders should be determined taking the root and comparison results into account
		        else{
			    hval=root;	
	         	}
	
	       	 if(bisectioncounter==15){           // for limiting loop
			 cout<<" The number of iterations for bisection method exceeded the threshold. " <<endl;
			 break;
	     	 }
            }
		}
		
		
		else if(signbit(f_lval)==signbit(f_hval)){
			cout<<"Warning! The signs of initial guesses are same. Not appropriate for bisection method."<<endl<<endl;
		} 
		
		else{
			
		bisectioncounter=0;
		
	    	while(hval-lval>tol){                // checking the convergence of boundary values of method
			root=lval+(hval-lval)/2;            // repeatedly halving the initial bracket until solution has been isolated
	    	f_hval=funct(array,argc,hval);     // getting f values of inital guesses and root;
		    f_lval=funct(array,argc,lval);
		    f_root=funct(array,argc,root);
		 
		    bisectioncounter++;
	        	if(signbit(f_lval)==signbit(f_root)){  
			    lval=root;
	        	}
	        	else{
			    hval=root;	
	         	}
		
	       	if(bisectioncounter==15){
			cout<<"The number of iterations for bisection method exceeded the threshold. " <<endl;
			break;
		    }
            }
}
    if(bisectioncounter!=0){
	cout<<"Root of bisection method: "<<root<<endl;
	cout<<"The number of iterations for bisection algorithm: "<<bisectioncounter<<endl<<endl;
}
}


void all::secant(float*array,float hval,float lval,float tol,float argc){
    float approot;    
	float check=-1;           // control variable for ending the loop
	int secantcounter=0;
   
	int flag=0;    // control variable for printing the results
	
	while(check==-1){	  
	float f_approot;
	secantcounter++;

	if(funct(array,argc,hval)-funct(array,argc,lval)==0){
		cout<<"Error! During the calculations in secant method one denominator became 0."<<endl<<endl;
		flag=1;        // result will have problems, could be NaN, shouldnt be printed, so control variable is changed
		break;
	}
	
	approot=hval-(funct(array,argc,hval)*(hval-lval)/(funct(array,argc,hval)-funct(array,argc,lval)));   // secant algorithm
	
    f_approot=funct(array,argc,approot);

	
	if(fabs(approot-hval)<tol){    // hval preserves the previous result of approot
    	check=0;                   // previous and new results are compared to see whether convergence is occured
	}
	
	lval=hval;      // determining the new border values
	hval=approot; 

	if(secantcounter==15){         // limiting the loop
		cout<<"The number of iterations for secant method exceeded the threshold. " <<endl;
		break;
	}
	}
	
	if(flag!=1){    // determining to print out the values or not via control variable
    cout<<"Root of secant method: "<<approot<<endl;
    cout<<"The number of iterations for secant algorithm: "<<secantcounter<<endl<<endl; }
}


void all::hybrid(float*array,float hval,float lval,float tol,float argc){
	  
	    int counter=0;
	    float f_lval,f_hval,f_root,root;  
	    
	    // bisection algoritm first , same with the above method except the loop limting counter
	    
	    f_hval=funct(array,argc,hval);     // getting f values of inital guesses and root;
	   	f_lval=funct(array,argc,lval);
			
	    if(f_lval==0 && f_hval==0){
				cout<<"Initial guesses are roots themselves. "<<endl<<endl;	//	break;
			}
			
		else if(f_lval==0 || f_hval==0 ){
				cout<<"One of the initial guesses is root itself. "<<endl<<endl;
			
		        root=lval+(hval-lval)/2;  
		
		    if(f_hval==0){
			hval=root;
	     	}
		    else if(f_lval==0){
			lval=root;
		    }
		
		    counter++;

	     	while(hval-lval>tol){    
		   
		    root=lval+(hval-lval)/2;  
	         counter++;
	    
	         	if(signbit(funct(array,argc,lval))==signbit(funct(array,argc,root))){  
			    lval=root;
	        	}
		        else{
			    hval=root;	
	         	}
	
	       	 if(counter==2){    // just 2 iteration with bisection method is desired
			 break;
	     	 }
            }
		}
		
		else if(signbit(f_lval)==signbit(f_hval)){
			cout<<"Warning! The signs of initial guesses are same. Not appropriate for hybrid method."<<endl<<endl;
		} 
		
		else{
			
		counter=0;
		
	    	while(hval-lval>tol){                // checking the convergence of boundary values of method
			root=lval+(hval-lval)/2;            // repeatedly halving the initial bracket until solution has been isolated
	    	f_hval=funct(array,argc,hval);     // getting f values of inital guesses and root;
		    f_lval=funct(array,argc,lval);
		    f_root=funct(array,argc,root);
		 
		    counter++;
	        	if(signbit(f_lval)==signbit(f_root)){  
			    lval=root;
	        	}
	        	else{
			    hval=root;	
	         	}
		
	       	if(counter==2){
			break;
		    }
            }
}
    // hval and lval already determined by bisection method
    
    float approot;    
    float check=-1;
    
	if(counter!=0){   // if bisection method isnt applicable, other part of the hybrid method(secant method) should also be terminated. Counter must be 2 if bisection method is used.
	while(check==-1){	  // rest of the secant algorithm is same with the method written above
	float f_approot;
	counter++;
	
	approot=hval-(funct(array,argc,hval)*(hval-lval)/(funct(array,argc,hval)-funct(array,argc,lval)));

    f_approot=funct(array,argc,approot); 

    if(fabs(approot-hval)<tol){
    	check=0;
	}
     
	lval=hval;
	hval=approot;
 
	if(counter==15){
		cout<<"The number of iterations for hybrid method exceeded the threshold. " <<endl;
		break;
	}
	   
}}

 
if(counter!=0){      // if bisection method isnt applicable, there would be no result the print out
cout<<"Root of hybrid method: "<<approot<<endl;
cout<<"The number of iterations in hybrid algorithm: "<<counter<<endl<<endl;
 }
}


int main( int argc, char *argv[] )  {                              
	float lowval,highval,tolerance,root;
	int coefficient_num=argc-4;   // number of coefficients of the polynomial
	lowval=atof(argv[1]);     // converting string variables to floats
	highval=atof(argv[2]);
	tolerance=atof(argv[3]);
    
    float*newarray;                 // dynamic allocation of an array to hold command line arguments(just coefficients of the polynomial)          
    newarray=new float[argc-4];
    
    for(int i=4;i<argc;i++){        // holding the coefficient of polynomial in array
    	newarray[i-4]=atof(argv[i]);
	}
	
	all methods;  // creating variable of class all
		methods.bisection(newarray,highval,lowval,tolerance,argc);
		methods.secant(newarray,highval,lowval,tolerance,argc);
        methods.hybrid(newarray,highval,lowval,tolerance,argc);
        
   delete newarray;     

   return 0;
}
