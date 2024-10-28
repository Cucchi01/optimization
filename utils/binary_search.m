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

function [pos] = binary_search(vector, val)
    % OUTPUT: the new value should be inserted at position pos+1
    L=1;
    n = length(vector);
    R = n;
    while L<=R
        m = floor((L+R)/2);
        if vector(m)==val
            pos = m;
            return;
        end
        if L==R
            break;
        end

        if vector(m) < val
            L = m+1;
        elseif vector(m) > val
            R = m-1;
        end
    end
    pos = L;
    %check if the new val should be after the last element
    if vector(pos)<val
        pos=pos+1;
    end  
end
