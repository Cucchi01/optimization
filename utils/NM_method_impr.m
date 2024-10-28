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

function [xmin,fmin,k,xseq] = NM_method_impr(f,x0,kmax,rho,chi,gamma,sigma,tol)
    %Function to apply to each xi at step k  in the shrinking phase
    fshrinking = @(x1,xi) x1 + sigma*(xi-x1) ;
    
    %Variable initialization
    k = 1;
    xk = x0;
    [m,n] = size(x0);
    % cr = 1;
    flag_shr = true;
    % xseq = zeros(m, 1);
    xseq = zeros(m,n,kmax);
    xseq(:,:,1) = x0;
    while k<kmax
        %SORTING PHASE: Calculate f(xi) for each xi in the symplex at step k, then
        %reorder the x in the symplex depending on the f(xi) value
        if flag_shr
            % a complete new sort should be done only at the start and after a shrinkage
            fxks = f(xk);
            [fxks,sortedInd] = sort(fxks);
            xk = xk(:,sortedInd);
            flag_shr=false;
        else
            fnew = f(xk(:,end));
            % the new position is found with a binary search
            % pos = binary_search(f(xk(:, 1:end-1)), fnew);
            % pos = find(f(xk(:, 1:end-1)) <= fnew, 1, 'last');
            % pos = position_insertion_sort(f(xk(:, 1:end-1)), fnew);
            % 
            % if pos <= 1
            %     xk = [xk(:, end), xk(:, 1:end-1)];
            %     fval = [fnew, fval(1:end-1)];
            % else%if not(pos >= n) % the check could be avoided because pos is equal to "num_points" only in the case of the Shrinkage
            %     xk = [xk(:, 1:pos-1), xk(:,end), xk(:, pos:end-1)];
            %     fval = [fval(1:pos-1), fval(end), fval(pos:end-1)];
            % end
            [xk,fxks ] = insertion_sort(xk, fxks(1:end), fnew);
        end
        change_rate = sqrt(sum((fxks-mean(fxks)).^2)/(n+1));
        if change_rate <= tol
            break
        end
        
        %REFLECTION PHASE  
        xb = sum(xk(:,1:end-1),2)/(n-1); %BARICENTER OF best "num_points"-1 points
        xr = xb + rho*(xb-xk(:,end));
        if(f(xr)<fxks(end-1) && f(xr)>= fxks(1))
            xk(:,end) = xr;
        elseif(f(xr)<fxks(1))
            %EXPANSION PHASE
            xe = xb + chi*(xr-xb);
            if(f(xe)<f(xr))
                xk(:,end) = xe;
                %xk(:,end-1) = xr;%scommentando questa linea si prende come "num_points"-1
                %esimo punto xr 
            else
                xk(:,end) = xr;
            end
        else
            %CONTRACTION PHASE
            if(f(xr)<fxks(end))
                xc = xb - gamma*(xb-xr);
                %xk(:,end) = xr; Commentando questa riga l'algoritmo esce come
                %lo inerpretavate voi, se la scommentate potete testare come
                %viene prendendo xr come xi in "num_points"
            else
                xc = xb - gamma*(xb-xk(:,end));
            end
            %xc = xb - gamma*(xb-xk(:,end));
            if(f(xc)>= fxks(end))
                %SHRINKAGE PHASE
                flag_shr=true;
                fshrSpecific = @(x)fshrinking(xk(:,1),x);% At each iteration the minimum changes and so it does the first parameter of the function                                                        
                xk(:,2:end) = fshrSpecific(xk(:,2:end));% It was possible to do all in one line without defining another function handler, i did it just to have a more readable code
            else
                xk(:,end) = xc;
            end
        end
        k = k+1;
        %cr = norm(xk-xkold);%change rate at each step just checks how the
        %symplex has changed over the iterations
        % fks = f(xk);
        % cr = sqrt(sum((fks-mean(fks)).^2)/(n+1));%si ferma qualche iterazione prima rispetto al cr di sopra
        xseq(:,:,k) = xk;
    end
    xseq = xseq(:,:,1:k);
    xmin = xk(:,1);
    fmin = f(xmin);
end