% Copyright 2024 Andrea Cucchietti, Davide Elio Stefano Demicheli, Shakti Singh Rathore
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%% Nelder-Mead on problem 3
clc

SEED = 123456;
rng(SEED);

n = 10;
delta = 0.3;
x0 = ones(n, 1);
vec_x0 = [ones(n, 1), 3*ones(n, 1), zeros(n, 1), rand(n, 1), linspace(-0.5, 0.5, n)']; 

banded_trig_grad = @(x, n) [(1:n-1)'.*sin(x(1:end-1, 1))+2*cos(x(1:end-1, 1)); n*sin(x(end, 1))-(n-1)*cos(x(end, 1))];
banded_trig_hess = @(x, n) diag([(1:n-1)'.*cos(x(1:end-1, 1))-2*sin(x(1:end-1, 1)); n*cos(x(end, :)) + (n-1)*sin(x(end,1))]);

f = @(x) banded_trig_general(x, n);
gradf = @(x) banded_trig_grad(x, n);
Hessf = @(x) banded_trig_hess(x, n);

tol = 1e-8;
kmax=100000;

%tuning
rho = 1;
chi = 1.5; %1.25
gamma = 0.9; % 0.8
sigma = 0.5;

% standard
%rho = 1;
%chi = 2;
%gamma = 0.5;
%sigma = 0.5;

for index_starting_point=1:size(vec_x0,2)

    X0 = repmat(vec_x0(:,index_starting_point), [1, n+1]) + [zeros(n, 1), delta* diag(ones(n,1))];
    
    disp("STARTING POINT:");
    disp(index_starting_point);
    disp("Nelder-Mead:")
    tic
    [xmin,fmin,k,xseq] = NM_method(f,X0,kmax,rho,chi,gamma,sigma,tol);
    toc
    disp("First 5 coordinates of the solution:")
    disp([mat2str(xmin(1:5))]);
    disp("Value of the solution:")
    disp(num2str(fmin));
    %disp(['xk: ', mat2str(xmin(1:5))]);
    %disp(['fk: ', num2str(fmin)]);
    %disp(['gradfk_norm: ', num2str(norm(gradf(xmin)), 2)]);
    %disp(['k: ', num2str(k)]);
    if k < kmax
        disp(["Converges in",num2str(k),"steps"])
    else
        disp("Does not converge")
    end
    
    
    disp("Nelder-Mead improved:")
    tic
    [xmin,fmin,k,xseq] = NM_method_impr(f,X0,kmax,rho,chi,gamma,sigma,tolgrad);
    toc
    disp("First 5 coordinates of the solution:")
    disp([mat2str(xmin(1:5))]);
    disp("Value of the solution:")
    disp(num2str(fmin));
    %disp(['gradfk_norm: ', num2str(norm(gradf(xmin)), 2)]);
    %disp(['k: ', num2str(k)]);
    if k < kmax
        disp(["Converges in",num2str(k),"steps"])
    else
        disp("Does not converge.")
    end

end
