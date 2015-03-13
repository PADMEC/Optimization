function [R,C,dR,D]=wscon(x,funame,coname,t,fm,l,a,n,par,varargin)

C=[];D=[];
global GradOutput_on
if GradOutput_on
    [R,~,dR]=feval(coname,x,varargin{:});
else
    R=feval(coname,x,varargin{:});
end