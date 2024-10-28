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

function [p] = compute_conv_rate(xseq, xsol)
    [~, k] = size(xseq);
    % k-2 because the last point in xseq is the solution
    y = zeros(1, k-2);
    x = zeros(1, k-2);

    for i=1:k-2
        y(1,i) = log(norm(xseq(:, i+1)- xsol)); %log||e_k+1||
        x(1,i) = log(norm(xseq(:, i) - xsol)); %log||e_k||
    end    

    p = x(:)\y(:);
end

