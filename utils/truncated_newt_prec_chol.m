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

function [xk, fk, gradfk_norm, k, xseq, btseq, internal_jseq, flag_converged, ...
    flag_armijo_always_satisfied, flag_enough_internal_steps] = truncated_newt_prec_chol(x0, f, gradf, Hessf, fterm, kmax, internal_jmax, ...
    tolgrad, c1, rho, btmax)

    flag_converged = false;
    flag_armijo_always_satisfied = true;
    flag_enough_internal_steps = true;
    
    farmijo = @(fk, alpha, gradfk, pk) ...
        fk + c1 * alpha * gradfk' * pk;
    
    dim_x0 = length(x0);
    xseq = zeros(dim_x0, kmax);
    btseq = zeros(1, kmax);
    internal_jseq = zeros(1, kmax);
    k=0;
    xk = x0;

    fk = f(xk);
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk, 2);

    while k<kmax && gradfk_norm >= tolgrad 
        zj = 0;
        rj = -gradfk;
        Bk = Hessf(xk);
        try 
            % M = L*L';
            flag_decomp = true;
            L = ichol(sparse(Bk));
        catch ME
            flag_decomp = false;
            L = eye(dim_x0);
        end
        if flag_decomp
            y0 = (L')\(L\rj);
        else
            y0 = rj;
        end
        yj = y0;
        %k is not useful for the forcing terms we have defined, but in
        %general it can be used
        tol = fterm(k, -rj) * gradfk_norm;
        dj = yj;
        j = 0;
        while j<internal_jmax
            if dj'* Bk *dj<=0
                if j ==0
                    % at the first step rj is -gradfk
                    pk = rj;
                    break;
                else
                    pk = zj;
                    break;
                end
            end
            w = Bk * dj;
            alpha = rj'*yj/(dj'*w);
            zj = zj + alpha * dj;
            rnew = rj - alpha * w;
            j = j+1;
            if norm(rnew, 2)<tol
                pk = zj;
                break;                
            end
            
            if flag_decomp
                ynew = (L')\(L\rnew);
            else
                ynew = rnew;
            end
            beta = rnew'*ynew/(rj'*yj);
            dj = ynew + beta *dj;
            rj = rnew;
            yj = ynew;
        end
        if  j>=internal_jmax && norm(rnew, 2)>=tol
            pk = zj;
            flag_enough_internal_steps = false;
        end

        alpha = 1;

        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        bt = 0;
        while bt < btmax && fnew >= farmijo(fk, alpha, gradfk, pk)
            alpha = rho * alpha;
            xnew = xk + alpha * pk;
            fnew = f(xnew);
            
            bt = bt + 1;            
        end

        if fnew >= farmijo(fk, alpha, gradfk, pk)
            flag_armijo_always_satisfied=false;
        end

        xk = xk + alpha * pk;

        fk = f(xk);
        gradfk = gradf(xk);
        gradfk_norm = norm(gradfk, 2);
        k=k+1;
        xseq(:, k) = xk;
        btseq(1, k) = bt;
        internal_jseq(1, k) = j;
    end
    
    if gradfk_norm < tolgrad
        flag_converged=true;
    end

    xseq = xseq(:, 1:k);
    btseq = btseq(:, 1:k);
    internal_jseq = internal_jseq(:, 1:k);

end

