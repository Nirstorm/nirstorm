function [idx, dist] = nst_knnsearch(X,Y,varargin)

% For compatibility
try
    [idx, dist] = knnsearch(X,Y,varargin{:});
    pp
catch
    % Alternate implementation from file exchange
    % https://www.mathworks.com/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit
    %
    % Licence:
    %     Copyright (c) 2009, Yi Cao
    %     All rights reserved.
    %
    %     Redistribution and use in source and binary forms, with or without
    %     modification, are permitted provided that the following conditions are met:
    %
    %     * Redistributions of source code must retain the above copyright notice, this
    %       list of conditions and the following disclaimer.
    %
    %     * Redistributions in binary form must reproduce the above copyright notice,
    %       this list of conditions and the following disclaimer in the documentation
    %       and/or other materials provided with the distribution
    %     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    %     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    %     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    %     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
    %     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    %     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    %     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    %     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    %     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    %     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    %
    % By Yi Cao at Cranfield University on 25 March 2008
    
    % Check inputs
    all_varagin = {Y X varargin{:}};
    [Q,R,K,fident] = parseinputs(all_varagin{:});
    
    % Check outputs
    nargoutchk(0,2);
    
    % C2 = sum(C.*C,2)';
    [N,M] = size(Q);
    L=size(R,1);
    idx = zeros(N,K);
    D = idx;
    
    if K==1
        % Loop for each query point
        for k=1:N
            d=zeros(L,1);
            for t=1:M
                d=d+(R(:,t)-Q(k,t)).^2;
            end
            if fident
                d(k)=inf;
            end
            [D(k),idx(k)]=min(d);
        end
    else
        for k=1:N
            d=zeros(L,1);
            for t=1:M
                d=d+(R(:,t)-Q(k,t)).^2;
            end
            if fident
                d(k)=inf;
            end
            [s,t]=sort(d);
            idx(k,:)=t(1:K);
            D(k,:)=s(1:K);
        end
    end
    
    if nargout>1
        dist=sqrt(D);
    end

end

end


function [Q,R,K,fident] = parseinputs(varargin)
% Copyright (c) 2009, Yi Cao
% See Licence in main function 

% Check input and output
narginchk(1,3);

Q=varargin{1};

if nargin<2
    R=Q;
    fident = true;
else
    fident = false;
    R=varargin{2};
end

if isempty(R)
    fident = true;
    R=Q;
end

if ~fident
    fident = isequal(Q,R);
end

if nargin<3
    K=1;
else
    K=varargin{3};
end

end