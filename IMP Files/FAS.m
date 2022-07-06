%% FAS Alg.
% Inputs:
% Psi - Approximation for SF                 W - Approximation for Vort.
% V - Approximation for azi. vel.            Sp - RHS of SF eqn                         
% Sw - RHS of Vort. eqn                      Sv - RHS of V eqn                          
% BC - Bound. Cond. function                 h - size of grid
% Re - Reynolds Number                       
% a - under-relaxation parameter             tol - tolerance
% Outputs:
% Psi - Solution for SF at grid size h      
% W - Solution for Vort. at grid size h
% V - Solution for azi. vel. at grid size h  
function [Psi,W,V] = FAS(Psi,W,V,Sp,Sw,Sv,BC,h,Re,a,tol)
    %% First Residual
    [Pres,Wres,Vres] = residual(Psi,W,V,Sp,Sw,Sv,h,Re); 
    r = max([abs(Pres);abs(Wres);abs(Vres)],[],'all');
    str = fprintf('residual = %f, %.2f%% Complete',[NaN,0]);
    
    %% Iterate until desired accuarcy
    while r>tol
        fprintf(repmat('\b',1,str))
        str = fprintf('residual = %f, %.2f%% Complete',[r,tol/r*100]);
        
        % V-Cycle
        [Psi,W,V] = Vcycle(Psi,W,V,Sp,Sw,Sv,BC,h,Re,a);

        % Residual 
        [Pres,Wres,Vres] = residual(Psi,W,V,Sp,Sw,Sv,h,Re); 
        r = max([abs(Pres);abs(Wres);abs(Vres)],[],'all');
    end
    fprintf(repmat('\b',1,str))
end