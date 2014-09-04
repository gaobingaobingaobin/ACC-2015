function [y, x, u] = trajectories(thetas, Ad, Bd, u_fun, x0, TR, N)
    %FUNCTION TRAJECTORIES computes trajectories of the model 
    %   Not meant to be public, only called from within toolbox functions 
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 
    
    n = size(Ad, 1); 
    x = zeros(n, N); 

    % compute input and state trajectories 
    u = u_fun(TR*(1:N));
    x(:, 1) = diag(cos(thetas(1,1:n))) * x0; 
    for k=1:N-1
        x(:, k+1) = Ad * diag(cos(thetas(k+1,1:n))) * x(:, k) + Bd * u(k); 
    end
    
    % concatenate state and input trajectories 
    y = sin(thetas').*[x; u];

end




