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

clear;
%clc;
close all;

load('forcing_term.mat');

SEED = 123456;
rng(SEED);
FILENAME_OUTPUT_EXCEL= 'results_point_3_P3.5.xlsx';

% set this based on nm
n = 1e1;
vec_n = [n, 1e2];

% standard
%rho = 0.5; 
rho = 0.8;
btmax = floor(log(1e-10)/log(rho));
c1 = 1e-4;
kmax = 100;
tolgrad = 1e-8;
%fterm = fterms_quad;

fterm_vec = cell(2,1);
fterm_vec{1} = fterms_suplin;
fterm_vec{2} = fterms_quad;
% fterm_vet{3} = fterms_lin;

% for the test on rho and c1
%c1_vec = [1e-4];
%[temp1, temp2] = meshgrid(rho_vec, c1_vec);
%rho_c1_cart_prod = [temp1(:), temp2(:)];

jmax = 100;
%one row of results per num n tested * num forcing terms tried * num starting
%point * num method
%each row has n, index_forcing_term, index_starting_point, index_method, if
%converged, armijo satisfied, bjmax never reached, k, conv_rate, dist_min, median btseq, mean 
% btseq, max btseq, median jseq, mean jseq, max jseq, time, fk, gradfk,
% the first 6 coeff of the final xk
results = zeros(length(vec_n)* length(fterm_vec)*5*2, 24);

row_res = 1;

disp('Truncated newton method all the alternatives: ')
for it_n=1:length(vec_n)
    disp("****************************")
    % rho = rho_c1_cart_prod(it_n, 1);
    % btmax = floor(log(1e-10)/log(rho));
    % c1 = rho_c1_cart_prod(it_n, 2);
    % disp(['RHO: ', num2str(rho)]);
    % disp(['SO WE HAVE BTMAX: ', num2str(btmax)]);
    % disp(['C1: ', num2str(c1)]);

    n = vec_n(it_n);

    % function handles of the function, its gradient and hessian
    banded_trig = @(x, n) sum((1:n)) - (1:n)*(cos(x(1:end,1))) + sum(2*sin(x(1:end-1, 1))) -(n-1)*sin(x(end, 1));
    banded_trig_grad = @(x, n) [(1:n-1)'.*sin(x(1:end-1, 1))+2*cos(x(1:end-1, 1)); n*sin(x(end, 1))-(n-1)*cos(x(end, 1))];
    banded_trig_hess = @(x, n) diag([(1:n-1)'.*cos(x(1:end-1, 1))-2*sin(x(1:end-1, 1)); n*cos(x(end, :)) + (n-1)*sin(x(end,1))]);
    
    f = @(x) banded_trig(x, n);
    gradf = @(x) banded_trig_grad(x, n);
    Hessf = @(x) banded_trig_hess(x, n);

    disp(['DIMENSION: ', num2str(n)]);
    % starting points
    vec_x0 = [ones(n, 1), 3*ones(n, 1), zeros(n, 1), rand(n, 1), linspace(-0.5,0.5,n)']; 

    for it_fterm = 1:length(fterm_vec)
        disp("---------------------------------------------")
        fterm = fterm_vec{it_fterm};
        disp("FORCING TERM: ")
        disp(func2str(fterm))
        
        [~, num_startin_x0] = size(vec_x0);
        for it_point = 1: num_startin_x0
            disp("+++++++++++++++++++++++++++++++++++++++++++++++")
            x0 = vec_x0(:, it_point);
            disp(["STARTING POINT INDEX AND SOME COEFF:", num2str(it_point)]);
            disp(x0(1:6));
            
            %row_res = row_res+1;
    
            disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            try 
                disp('TRUNCATED NEWTON METHOD:')
                tic
                [xk, fk, gradfk_norm, k, xseq, btseq, jseq, flag_converged,...
                    flag_always_armijo_satisfied, flag_enough_internal_steps] = truncated_newt(x0, f, gradf, Hessf, fterm, kmax, jmax, ...
                tolgrad, c1, rho, btmax);
                time = toc;
                disp("Time needed:")
                disp(time)
                disp(["Converged: ", mat2str(flag_converged)]);
                disp(['fk: ', num2str(fk),';  gradfk_norm: ', num2str(gradfk_norm), ';  xk: ', mat2str(xk)]);
                disp(['k: ', num2str(k)]);
                disp(["Armijo always satisfied: ", mat2str(flag_always_armijo_satisfied)]);
                disp('btseq: ');
                disp(size(btseq));
                disp(["Median btseq: ", num2str(median(btseq)), ", Mean btseq: ", ...
                    num2str(mean(btseq)), ";  Max btseq: ", num2str(max(btseq))]);
                % disp(btseq);
                disp('xseq: ');
                disp(size(xseq));  
                disp(["Never reached the maximum number of internal steps: ", mat2str(flag_enough_internal_steps)]);
                disp('jseq: ');
                disp(size(jseq));
                disp(["Median jseq: ", num2str(median(jseq)), ";  Mean jseq: ", num2str(mean(jseq)), ";  Max jseq: ", num2str(max(jseq))]);
                
                if flag_converged == true
                    sol = xseq(:, end); 
                    rate_conv = compute_conv_rate(xseq, sol);
                    disp(["Conv. rate", num2str(rate_conv)]);
                else 
                    rate_conv = 0;
                end
            
                %index of "normal" Truncated method is 1
                results(row_res, :) = [n, it_fterm, it_point, 1, flag_converged,...
                    flag_always_armijo_satisfied, flag_enough_internal_steps,k, rate_conv,...
                    median(btseq), mean(btseq), max(btseq), median(jseq), mean(jseq), max(jseq), time, fk, gradfk_norm, xk(1:6)'];
                        
            catch ME
                disp(ME.identifier);
                if ME.identifier == "MATLAB:array:SizeLimitExceeded"
                    results(row_res, 5) = "SIZELIMITERROR";
                    disp("SIZE LIMIT EXCEEDED");
                else
                    results(row_res, 5) = "ERROR";
                    rethrow(ME);
                end
            end
            row_res = row_res+1;
            
            disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            try
                disp('TRUNCATED NEWTON METHOD PRECONDITIONING:')
                tic
                [xk, fk, gradfk_norm, k, xseq, btseq, jseq, flag_converged,...
                    flag_always_armijo_satisfied, flag_enough_internal_steps] = truncated_newt_prec(x0, f, gradf, Hessf, fterm, kmax, jmax, ...
                tolgrad, c1, rho, btmax);
                time = toc;
                disp("Time needed:")
                disp(time)

                disp(["Converged: ", mat2str(flag_converged)]);
                disp(['fk: ', num2str(fk),';  gradfk_norm: ', num2str(gradfk_norm), ';  xk: ', mat2str(xk)]);
                disp(['k: ', num2str(k)]);
                disp(["Armijo always satisfied: ", mat2str(flag_always_armijo_satisfied)]);
                disp('btseq: ');
                disp(size(btseq));
                disp(["Median btseq: ", num2str(median(btseq)), ", Mean btseq: ", ...
                    num2str(mean(btseq)), ";  Max btseq: ", num2str(max(btseq))]);
                % disp(btseq);
                disp('xseq: ');
                disp(size(xseq));  
                disp(["Never reached the maximum number of internal steps: ", mat2str(flag_enough_internal_steps)]);
                disp('jseq: ');
                disp(size(jseq));
                disp(["Median jseq: ", num2str(median(jseq)), ";  Mean jseq: ", num2str(mean(jseq)), ";  Max jseq: ", num2str(max(jseq))]);
                
                sol = xseq(:, end);
                rate_conv = compute_conv_rate(xseq, sol);
                disp(["Conv. rate", num2str(rate_conv)])
            
                %index of preconditioned Truncated method is 2
                results(row_res, :) = [n, it_fterm, it_point, 2, flag_converged,...
                    flag_always_armijo_satisfied, flag_enough_internal_steps,k, rate_conv,...
                    median(btseq), mean(btseq), max(btseq), median(jseq), mean(jseq), max(jseq), time, fk, gradfk_norm, xk(1:6)'];
                %end
            
            catch ME
                disp(ME.identifier);
                if ME.identifier == "MATLAB:array:SizeLimitExceeded"
                    disp("SIZE LIMIT EXCEEDED");
                else
                    rethrow(ME);
                end
            end
            row_res = row_res+1;
        end
    end
end

writematrix(results,FILENAME_OUTPUT_EXCEL,'Sheet',1,'Range','A2:AC100');
disp("END");
