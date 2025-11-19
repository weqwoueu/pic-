SUBROUTINE  Retangular_reference_basis(x,y,basis_type,basis_index,derivative_degree_x,derivative_degree_y,r)

IMPLICIT NONE

REAL(8)                                               x, y !, kesai, eita

INTEGER                                               basis_type, basis_index, derivative_degree_x, derivative_degree_y

REAL(8)                                               r

IF (basis_type==1) THEN
    
    IF (derivative_degree_x==0 .AND. derivative_degree_y==0) THEN
        
        IF (basis_index==1) THEN
            r=(-x-y+x*y+1)/4
        ELSEIF (basis_index==2) THEN
            r=(x-y-x*y+1)/4
        ELSEIF (basis_index==3) THEN
            r=(x+y+x*y+1)/4
        ELSEIF (basis_index==4) THEN
            r=(-x+y-x*y+1)/4
        ENDIF
              
    ELSEIF (derivative_degree_x==1 .AND. derivative_degree_y==0) THEN
            
        IF (basis_index==1) THEN
            r=(-1+y)/4
        ELSEIF (basis_index==2) THEN
            r=(1-y)/4
        ELSEIF (basis_index==3) THEN
            r=(1+y)/4
        ELSEIF (basis_index==4) THEN
            r=(-1-y)/4
        ENDIF
       
    ELSEIF (derivative_degree_x==0 .AND. derivative_degree_y==1) THEN
            
        IF (basis_index==1) THEN
            r=(-1+x)/4
        ELSEIF (basis_index==2) THEN
            r=(-1-x)/4
        ELSEIF (basis_index==3) THEN
            r=(1+x)/4
        ELSEIF (basis_index==4) THEN
            r=(1-x)/4
        ENDIF
            
    ELSEIF (derivative_degree_x==1 .AND. derivative_degree_y==1) THEN
 
        IF (basis_index==1) THEN
            r=1/4+x-x
        ELSEIF (basis_index==2) THEN
            r=-1/4+x-x;
        ELSEIF (basis_index==3) THEN
            r=1/4+x-x;
        ELSEIF (basis_index==4) THEN
            r=-1/4+x-x;
        ENDIF
        
    ENDIF
    
ELSEIF (basis_type==2) THEN
    
    IF (derivative_degree_x==0 .AND. derivative_degree_y==0) THEN
        
        IF (basis_index==1) THEN
            r=(x*x*y*y-x*x*y-x*y*y+x*y)/4;
        ELSEIF (basis_index==2) THEN
            r=(x*x*y*y-x*x*y+x*y*y-x*y)/4;
        ELSEIF (basis_index==3) THEN
             r=(x*x*y*y+x*x*y+x*y*y+x*y)/4;
        ELSEIF (basis_index==4) THEN
             r=(x*x*y*y+x*x*y-x*y*y-x*y)/4;
        ELSEIF (basis_index==5) THEN
             r=(-x*x*y*y+x*x*y+y*y-y)/2;
        ELSEIF (basis_index==6) THEN
             r=(-x*x*y*y-x*y*y+x*x+x)/2;
        ELSEIF (basis_index==7) THEN
             r=(-x*x*y*y-x*x*y+y*y+y)/2;
        ELSEIF (basis_index==8) THEN
             r=(-x*x*y*y+x*y*y+x*x-x)/2;
        ELSEIF (basis_index==9) THEN
             r=x*x*y*y-x*x-y*y+1;
        ENDIF
             
    ELSEIF (derivative_degree_x==1 .AND. derivative_degree_y==0) THEN
            
        IF  (basis_index==1) THEN
             r=(2.*x*y*y-2*x*y-y*y+y)/4;
        ELSEIF (basis_index==2) THEN
             r=(2*x*y*y-2*x*y+y*y-y)/4;
        ELSEIF (basis_index==3) THEN
             r=(2*x*y*y+2*x*y+y*y+y)/4;
        ELSEIF (basis_index==4) THEN
             r=(2*x*y*y+2*x*y-y*y-y)/4;
        ELSEIF (basis_index==5) THEN
             r=(-2*x*y*y+2*x*y)/2;
        ELSEIF (basis_index==6) THEN
             r=(-2*x*y*y-y*y+2*x+1)/2;
        ELSEIF (basis_index==7) THEN
             r=(-2*x*y*y-2*x*y)/2;
        ELSEIF (basis_index==8) THEN
             r=(-2*x*y*y+y*y+2*x-1)/2;
        ELSEIF (basis_index==9) THEN
             r=2*x*y*y-2*x;
        ENDIF
                      
    ELSEIF (derivative_degree_x==0 .AND. derivative_degree_y==1) THEN
            
        IF (basis_index==1) THEN
             r=(2*x*x*y-x*x-2*x*y+x)/4;
        ELSEIF (basis_index==2) THEN
             r=(2*x*x*y-x*x+2*x*y-x)/4;
        ELSEIF (basis_index==3) THEN
             r=(2*x*x*y+x*x+2*x*y+x)/4;
        ELSEIF (basis_index==4) THEN
             r=(2*x*x*y+x*x-2*x*y-x)/4;
        ELSEIF (basis_index==5) THEN
             r=(-2*x*x*y+x*x+2*y-1)/2;
        ELSEIF (basis_index==6) THEN
             r=(-2*x*x*y-2*x*y)/2;
        ELSEIF (basis_index==7) THEN
             r=(-2*x*x*y-x*x+2*y+1)/2;
        ELSEIF (basis_index==8) THEN
             r=(-2*x*x*y+2*x*y)/2;
        ELSEIF (basis_index==9) THEN
             r=2*x*x*y-2*y;
        ENDIF
      
    ELSEIF (derivative_degree_x==1 .AND. derivative_degree_y==1) THEN  
        
        IF (basis_index==1) THEN
             r=(4*x*y-2*x-2*y+1)/4;
        ELSEIF (basis_index==2) THEN
             r=(4*x*y-2*x+2*y-1)/4;
        ELSEIF (basis_index==3) THEN
             r=(4*x*y+2*x+2*y+1)/4;
        ELSEIF (basis_index==4) THEN
             r=(4*x*y+2*x-2*y-1)/4;
        ELSEIF (basis_index==5) THEN
             r=(-4*x*y+2*x)/2;
        ELSEIF (basis_index==6) THEN
             r=(-4*x*y-2*y)/2;
        ELSEIF (basis_index==7) THEN
             r=(-4*x*y-2*x)/2;
        ELSEIF (basis_index==8) THEN
             r=(-4*x*y+2*y)/2;
        ELSEIF (basis_index==9) THEN
             r=4*x*y;
        ENDIF
      
    ENDIF
        
ENDIF

END