function [xm,Ftrue,Fquadra]=SliceF_Quadratic(func,m,Cov_m,CC,index)
%
% this function computes the function and its quadratic approximation around
% a given value of the function parameter along a 1D slice of the 
% multi-dimensional space
%
% Inputs:
%   func  :: true function
%   m     :: points where the slice will be centered.
%   Cov   :: covariance matrix (invert of Hessian)
%   index :: dimension on which the slice is performed
%

q=length(m);
xcenter=m(index);

%% cholesky of covariance....
%% 
st_x=sqrt(CC(index,index));

half_seg_size=(1*st_x);

half_nsamples=25;
step=half_seg_size/half_nsamples;

xm=[-half_seg_size:step:half_seg_size]'+m(index);

F_at_m=feval(func,m);

Hess=inv(Cov_m);

for i=1:length(xm)
    
    maux=m;
    maux(index)=xm(i);
    
    Ftrue(i)=feval(func,maux);
    
    Fquadra(i)=0.5*(m-maux)'*Hess*(m-maux)+F_at_m;
    
    
end
   % Fquadra=Fquadra;
    
