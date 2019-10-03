function [W,M,V,L,E] = EM_GM(X,k,ltol,maxiter,pflag,Init)  
% Inputs:  
% X(n,d) - input data, n=number of observations, d=dimension of variable  
% k - maximum number of Gaussian components allowed  
% ltol - percentage of the log likelihood difference between 2 iterations ([] for none)  
% maxiter - maximum number of iteration allowed ([] for none)  
% pflag - 1 for plotting GM for 1D or 2D cases only, 0 otherwise ([] for none)  
% Init - structure of initial W, M, V: Init.W, Init.M, Init.V ([] for none)  
%  
% Ouputs:  
% W(1,k) - estimated weights of GM  
% M(d,k) - estimated mean vectors of GM  
% V(d,d,k) - estimated covariance matrices of GM  
% L - log likelihood of estimates  
%   
 
% % % % Validate inputs %%%%  
if nargin <= 1,  
    disp('EM_GM must have at least 2 inputs: X,k!/n')  
    return  
elseif nargin == 2,  
    ltol = 0.1; maxiter = 1000; pflag = 0; Init = [];  
    err_X = Verify_X(X);  
    err_k = Verify_k(k);  
    if err_X | err_k, return; end  
elseif nargin == 3,  
    maxiter = 1000; pflag = 0; Init = [];  
    err_X = Verify_X(X);  
    err_k = Verify_k(k);  
    [ltol,err_ltol] = Verify_ltol(ltol);  
    if err_X | err_k | err_ltol, return; end  
elseif nargin == 4,  
    pflag = 0; Init = [];  
    err_X = Verify_X(X);  
    err_k = Verify_k(k);  
    [ltol,err_ltol] = Verify_ltol(ltol);  
    [maxiter,err_maxiter] = Verify_maxiter(maxiter);  
    if err_X | err_k | err_ltol | err_maxiter, return; end  
elseif nargin == 5,  
     Init = [];  
    err_X = Verify_X(X);  
    err_k = Verify_k(k);  
    [ltol,err_ltol] = Verify_ltol(ltol);  
    [maxiter,err_maxiter] = Verify_maxiter(maxiter);  
    [pflag,err_pflag] = Verify_pflag(pflag);  
    if err_X | err_k | err_ltol | err_maxiter | err_pflag, return; end  
elseif nargin == 6,  
    err_X = Verify_X(X);  
    err_k = Verify_k(k);  
    [ltol,err_ltol] = Verify_ltol(ltol);  
    [maxiter,err_maxiter] = Verify_maxiter(maxiter);  
    [pflag,err_pflag] = Verify_pflag(pflag);  
    [Init,err_Init]=Verify_Init(Init);  
    if err_X | err_k | err_ltol | err_maxiter | err_pflag | err_Init, return; end  
else  
    disp('EM_GM must have 2 to 6 inputs!');  
    return  
end  
 
% % % % Initialize W, M, V,L %%%%  
t = cputime;  
if isempty(Init),  
    [W,M,V] = Init_EM(X,k); L = 0;  
else  
    W = Init.W;  
    M = Init.M;  
    V = Init.V;  
end  
[E,Ln] = Expectation(X,k,W,M,V); % Initialize log likelihood and compute expectation  
Lo = 2*Ln;  
 
% % % % EM algorithm %%%%  
niter = 0;  
while (abs(100*(Ln-Lo)/Lo)>ltol) && (niter<=maxiter),  
%**************************************************************************  
% For the first loop, expectation is computed above (cf. Line 69). For the  
% next loops, it is computed in their preceding loops. Therefore, the  
% following line [Line 80] is omitted.  
%**************************************************************************  
% E = Expectation(X,k,W,M,V); % E-step  
 
    [W,M,V] = Maximization(X,k,E); % M-step  
    Lo = Ln;  
    [E,Ln] = Expectation(X,k,W,M,V); % Log-likelihood and expectation computation  
    niter = niter + 1;  
    LnValues(niter) = Ln;  
    %disp(['Iteration # = ' num2str(niter) ', Ln = ' num2str(Ln) ', Stopping Criterion = ' num2str(abs(100*(Ln-Lo)/Lo)) ', Tolerance = ' num2str(ltol)])  
end  
L = Ln;  
% % % % Plot 1D or 2D %%%%  
if pflag==1,  
    [n,d] = size(X);  
    if d>2,  
        disp('Can only plot 1 or 2 dimensional applications!/n');  
    else  
        Plot_GM(X,k,W,M,V);  
    end  
    elapsed_time = sprintf('CPU time used for EM_GM: %5.2fs',cputime-t);  
    %disp(elapsed_time);  
    %disp(sprintf('Number of iterations: %d',niter-1));  
end  
%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of EM_GM %%%%  
%%%%%%%%%%%%%%%%%%%%%%  
 
function [E,L] = Expectation(X,k,W,M,V)  
[n,d] = size(X);  
a = (2*pi)^(0.5*d);  
S = zeros(1,k);  
iV = zeros(d,d,k);  
for j=1:k,  
    if V(:,:,j)==zeros(d,d), V(:,:,j)=ones(d,d)*eps; end  
    S(j) = sqrt(det(V(:,:,j)));  
    iV(:,:,j) = inv(V(:,:,j));  
end  
%*************************************************************************
%*************************************************************************  
E = zeros(n,k);  
chunkSize = 1000;  
howManyFullChunks = fix(n/chunkSize);  
numberOfRemainingVectors = n - howManyFullChunks*chunkSize;  
pdfValueVector = zeros(n,1);  
for j = 1:k  
    for chunkCounter = 1:howManyFullChunks  
        modificationRange = (chunkCounter-1)*chunkSize+1:chunkCounter*chunkSize;  
        meanSubtractedChunk = X(modificationRange,:)'-repmat(M(:,j),1,chunkSize);  
        pdfValueVectorToExp = diag(-0.5*meanSubtractedChunk'*iV(:,:,j)*meanSubtractedChunk);  
        pdfValueVector(modificationRange) = exp(pdfValueVectorToExp)/(a*S(j));  
    end  
    modificationRange = howManyFullChunks*chunkSize+1:n;  
    meanSubtractedChunk = X(modificationRange,:)'-repmat(M(:,j),1,length(modificationRange));  
    pdfValueVectorToExp = diag(-0.5*meanSubtractedChunk'*iV(:,:,j)*meanSubtractedChunk);  
    pdfValueVector(modificationRange) = exp(pdfValueVectorToExp)/(a*S(j));  
    E(:,j) = W(j)*pdfValueVector';  
end  
sumInColumnDirection = sum(E,2);  
divisor = repmat(sumInColumnDirection,1,k);  
L = sum(log(sum(E,2)));  
E = E./divisor;  
%**************************************************************************  
%*************************** MODIFICATION ENDS HERE ***********************  
%**************************************************************************  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Expectation %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
function [W,M,V] = Maximization(X,k,E)  
[n,d] = size(X);  
V = zeros(d,d,k);  
W = sum(E,1);  
M = X'*E;  
M = M./repmat(W,d,1);  

for i=1:k,  
    for j=1:n,  
        dXM = X(j,:)'-M(:,i);  
        V(:,:,i) = V(:,:,i) + E(j,i)*(dXM*dXM');
    end  
    V(:,:,i) = V(:,:,i)/W(i);  
end  
W = W/n;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Maximization %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
function err_X = Verify_X(X)  
err_X = 1;  
[n,d] = size(X);  
if n<d,  
    disp('Input data must be n x d!/n');  
    return  
end  
err_X = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Verify_X %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%  
 
function err_k = Verify_k(k)  
err_k = 1;  
if ~isnumeric(k) | ~isreal(k) | k<1,  
    disp('k must be a real integer >= 1!/n');  
    return  
end  
err_k = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Verify_k %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%  
 
function [ltol,err_ltol] = Verify_ltol(ltol)  
err_ltol = 1;  
if isempty(ltol),  
    ltol = 1e-2;  
elseif ~isreal(ltol) | ltol<=0,  
    disp('ltol must be a positive real number!');  
    return  
end  
err_ltol = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Verify_ltol %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
function [maxiter,err_maxiter] = Verify_maxiter(maxiter)  
err_maxiter = 1;  
if isempty(maxiter),  
    maxiter = 1000;  
elseif ~isreal(maxiter) | maxiter<=0,  
    disp('ltol must be a positive real number!');  
    return  
end  
err_maxiter = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Verify_maxiter %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
function [pflag,err_pflag] = Verify_pflag(pflag)  
err_pflag = 1;  
if isempty(pflag),  
    pflag = 0;  
elseif pflag~=0 & pflag~=1,  
    disp('Plot flag must be either 0 or 1!/n');  
    return  
end  
err_pflag = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Verify_pflag %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
function [Init,err_Init] = Verify_Init(Init)  
err_Init = 1;  
if isempty(Init),  
    % Do nothing;  
elseif isstruct(Init),  
    [Wd,Wk] = size(Init.W);  
    [Md,Mk] = size(Init.M);  
    [Vd1,Vd2,Vk] = size(Init.V);  
    if Wk~=Mk | Wk~=Vk | Mk~=Vk,  
        disp('k in Init.W(1,k), Init.M(d,k) and Init.V(d,d,k) must equal!/n')  
        return  
    end  
    if Md~=Vd1 | Md~=Vd2 | Vd1~=Vd2,  
        disp('d in Init.W(1,k), Init.M(d,k) and Init.V(d,d,k) must equal!/n')  
        return  
    end  
else  
    disp('Init must be a structure: W(1,k), M(d,k), V(d,d,k) or []!');  
    return  
end  
err_Init = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Verify_Init %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
function [W,M,V] = Init_EM(X,k)  
[n,d] = size(X);  
[Ci,C] = mykmeans(X,k,'Start','cluster', ...  
'Maxiter',100, ...  
'EmptyAction','drop', ...  
'Display','off'); % Ci(nx1) - cluster indeices; C(k,d) - cluster centroid (i.e. mean)  
while sum(isnan(C))>0,  
[Ci,C] = mykmeans(X,k,'Start','cluster', ...  
'Maxiter',100, ...  
'EmptyAction','drop', ...  
'Display','off');  
end  
M = C';  
Vp = repmat(struct('count',0,'X',zeros(n,d)),1,k);  
for i=1:n, % Separate cluster points  
Vp(Ci(i)).count = Vp(Ci(i)).count + 1;  
Vp(Ci(i)).X(Vp(Ci(i)).count,:) = X(i,:);  
end  
V = zeros(d,d,k);  
 
% Preallocating W  
W = zeros(1,k);  
 
for i=1:k,  
W(i) = Vp(i).count/n;  
V(:,:,i) = cov(Vp(i).X(1:Vp(i).count,:));  
end  
%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Init_EM %%%%  
%%%%%%%%%%%%%%%%%%%%%%%% 
 
function Plot_GM(X,k,W,M,V)  
[n,d] = size(X);  
if d>2,  
    disp('Can only plot 1 or 2 dimensional applications!/n');  
    return  
end  
S = zeros(d,k);  
R1 = zeros(d,k);  
R2 = zeros(d,k);  
for i=1:k, % Determine plot range as 4 x standard deviations  
    S(:,i) = sqrt(diag(V(:,:,i)));  
    R1(:,i) = M(:,i)-4*S(:,i);  
    R2(:,i) = M(:,i)+4*S(:,i);  
end  
Rmin = min(min(R1));  
Rmax = max(max(R2));  
R = [Rmin:0.001*(Rmax-Rmin):Rmax];  
% clf, hold on  
figure; hold on  
if d==1,  
    Q = zeros(size(R));  
    for i=1:k,  
        P = W(i)*normpdf(R,M(:,i),sqrt(V(:,:,i)));  
        Q = Q + P;  
        plot(R,P,'r-'); grid on,  
    end  
    plot(R,Q,'k-');  
    xlabel('X');  
    ylabel('Probability density');  
else % d==2  
% plot(X(:,1),X(:,2),'r.');  
% %     plot(X(1:200,1),X(1:200,2),'r.');  
% %     plot(X(201:400,1),X(201:400,2),'g.');  
% %     plot(X(401:600,1),X(401:600,2),'b.');  
    plot(X(:,1),X(:,2),'r.')
    for i=1:k,  
        Plot_Std_Ellipse(M(:,i),V(:,:,i));  
    end  
    xlabel('1^{st} dimension');  
    ylabel('2^{nd} dimension');  
    axis([Rmin Rmax Rmin Rmax])  
end  
title('Gaussian Mixture estimated by EM');  
%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Plot_GM %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%  
 
function Plot_Std_Ellipse(M,V)  
[Ev,D] = eig(V);  
d = length(M);  
if V(:,:)==zeros(d,d),
    disp('Warning Variance is null'), return
% %    V(:,:) = ones(d,d)*eps;  
end  
iV = inv(V);  
% Find the larger projection  
P = [1,0;0,0]; % X-axis projection operator  
P1 = P * 2*sqrt(D(1,1)) * Ev(:,1);  
P2 = P * 2*sqrt(D(2,2)) * Ev(:,2);  
if abs(P1(1)) >= abs(P2(1)),  
    Plen = P1(1);  
else  
    Plen = P2(1);  
end  
count = 1;  
step = 0.001*Plen;  
Contour1 = zeros(2001,2);  
Contour2 = zeros(2001,2);  
for x = -Plen:step:Plen,  
    a = iV(2,2);  
    b = x * (iV(1,2)+iV(2,1));  
    c = (x^2) * iV(1,1) - 1;  
    Root1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);  
    Root2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);  
    if isreal(Root1),  
        Contour1(count,:) = [x,Root1] + M';  
        Contour2(count,:) = [x,Root2] + M';  
        count = count + 1;  
    end  
end  
Contour1 = Contour1(1:count-1,:);  
Contour2 = [Contour1(1,:);Contour2(1:count-1,:);Contour1(count-1,:)];  
plot(M(1),M(2),'k+');  
plot(Contour1(:,1),Contour1(:,2),'k-');  
plot(Contour2(:,1),Contour2(:,2),'k-');  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% End of Plot_Std_Ellipse %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  




function [idx, C, sumD, D] = mykmeans(X, k, varargin)
%KMEANS K-means clustering.
%   IDX = KMEANS(X, K) partitions the points in the N-by-P data matrix
%   X into K clusters.  This partition minimizes the sum, over all
%   clusters, of the within-cluster sums of point-to-cluster-centroid
%   distances.  Rows of X correspond to points, columns correspond to
%   variables.  KMEANS returns an N-by-1 vector IDX containing the
%   cluster indices of each point.  By default, KMEANS uses squared
%   Euclidean distances.
%
%   KMEANS treats NaNs as missing data, and removes any rows of X that
%   contain NaNs. 
%
%   [IDX, C] = KMEANS(X, K) returns the K cluster centroid locations in
%   the K-by-P matrix C.
%
%   [IDX, C, SUMD] = KMEANS(X, K) returns the within-cluster sums of
%   point-to-centroid distances in the 1-by-K vector sumD.
%
%   [IDX, C, SUMD, D] = KMEANS(X, K) returns distances from each point
%   to every centroid in the N-by-K matrix D.
%
%   [ ... ] = KMEANS(..., 'PARAM1',val1, 'PARAM2',val2, ...) allows you to
%   specify optional parameter name/value pairs to control the iterative
%   algorithm used by KMEANS.  Parameters are:
%
%   'Distance' - Distance measure, in P-dimensional space, that KMEANS
%      should minimize with respect to.  Choices are:
%            {'sqEuclidean'} - Squared Euclidean distance
%             'cityblock'    - Sum of absolute differences, a.k.a. L1
%             'cosine'       - One minus the cosine of the included angle
%                              between points (treated as vectors)
%             'correlation'  - One minus the sample correlation between
%                              points (treated as sequences of values)
%             'Hamming'      - Percentage of bits that differ (only
%                              suitable for binary data)
%
%   'Start' - Method used to choose initial cluster centroid positions,
%      sometimes known as "seeds".  Choices are:
%                 {'sample'} - Select K observations from X at random
%                  'uniform' - Select K points uniformly at random from
%                              the range of X.  Not valid for Hamming distance.
%                  'cluster' - Perform preliminary clustering phase on
%                              random 10% subsample of X.  This preliminary
%                              phase is itself initialized using 'sample'.
%                  matrix    - A K-by-P matrix of starting locations.  In
%                              this case, you can pass in [] for K, and
%                              KMEANS infers K from the first dimension of
%                              the matrix.  You can also supply a 3D array,
%                              implying a value for 'Replicates'
%                              from the array's third dimension.
%
%   'Replicates' - Number of times to repeat the clustering, each with a
%      new set of initial centroids [ positive integer | {1}]
%
%   'Maxiter' - The maximum number of iterations [ positive integer | {100}]
%
%   'EmptyAction' - Action to take if a cluster loses all of its member
%      observations.  Choices are:
%               {'error'}    - Treat an empty cluster as an error
%                'drop'      - Remove any clusters that become empty, and
%                              set corresponding values in C and D to NaN.
%                'singleton' - Create a new cluster consisting of the one
%                              observation furthest from its centroid.
%
%   'Display' - Display level [ 'off' | {'notify'} | 'final' | 'iter' ]
%
%   Example:
%
%       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
%       [cidx, ctrs] = kmeans(X, 2, 'dist','city', 'rep',5, 'disp','final');
%       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%
%   See also LINKAGE, CLUSTERDATA, SILHOUETTE.

%   KMEANS uses a two-phase iterative algorithm to minimize the sum of
%   point-to-centroid distances, summed over all K clusters.  The first
%   phase uses what the literature often describes as "batch" updates,
%   where each iteration consists of reassigning points to their nearest
%   cluster centroid, all at once, followed by recalculation of cluster
%   centroids. This phase may be thought of as providing a fast but
%   potentially only approximate solution as a starting point for the
%   second phase.  The second phase uses what the literature often
%   describes as "on-line" updates, where points are individually
%   reassigned if doing so will reduce the sum of distances, and cluster
%   centroids are recomputed after each reassignment.  Each iteration
%   during this second phase consists of one pass though all the points.
%   KMEANS can converge to a local optimum, which in this case is a
%   partition of points in which moving any single point to a different
%   cluster increases the total sum of distances.  This problem can only be
%   solved by a clever (or lucky, or exhaustive) choice of starting points.
%
% References:
%
%   [1] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.
%   [2] Spath, H. (1985) Cluster Dissection and Analysis: Theory, FORTRAN
%       Programs, Examples, translated by J. Goldschmidt, Halsted Press,
%       New York, 226 pp.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.7 $  $Date: 2005/11/18 14:28:09 $

if nargin < 2
    error('stats:kmeans:TooFewInputs','At least two input arguments required.');
end

if any(isnan(X(:)))
    warning('stats:kmeans:MissingDataRemoved','Removing rows of X with missing data.');
    X = X(~any(isnan(X),2),:);
end

% n points in p dimensional space
[n, p] = size(X);
Xsort = []; Xord = [];

pnames = {   'distance'  'start' 'replicates' 'maxiter' 'emptyaction' 'display'};
dflts =  {'sqeuclidean' 'sample'          []       100        'error'  'notify'};
[eid,errmsg,distance,start,reps,maxit,emptyact,display] ...
                       = statgetargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('stats:kmeans:%s',eid),errmsg);
end

if ischar(distance)
    distNames = {'sqeuclidean','cityblock','cosine','correlation','hamming'};
    i = strmatch(lower(distance), distNames);
    if length(i) > 1
        error('stats:kmeans:AmbiguousDistance', ...
              'Ambiguous ''distance'' parameter value:  %s.', distance);
    elseif isempty(i)
        error('stats:kmeans:UnknownDistance', ...
              'Unknown ''distance'' parameter value:  %s.', distance);
    end
    distance = distNames{i};
    switch distance 
    case 'cityblock'
        [Xsort,Xord] = sort(X,1);
    case 'cosine'
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error('stats:kmeans:ZeroDataForCos', ...
                  ['Some points have small relative magnitudes, making them ', ...
                   'effectively zero.\nEither remove those points, or choose a ', ...
                   'distance other than ''cosine''.']);
        end
        X = X ./ Xnorm(:,ones(1,p));
    case 'correlation'
        X = X - repmat(mean(X,2),1,p);
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error('stats:kmeans:ConstantDataForCorr', ...
                  ['Some points have small relative standard deviations, making them ', ...
                   'effectively constant.\nEither remove those points, or choose a ', ...
                   'distance other than ''correlation''.']);
        end
        X = X ./ Xnorm(:,ones(1,p));
    case 'hamming'
        if ~all(ismember(X(:),[0 1]))
            error('stats:kmeans:NonbinaryDataForHamm', ...
                  'Non-binary data cannot be clustered using Hamming distance.');
        end
    end
else
    error('stats:kmeans:InvalidDistance', ...
          'The ''distance'' parameter value must be a string.');
end

if ischar(start)
    startNames = {'uniform','sample','cluster'};
    i = strmatch(lower(start), startNames);
    if length(i) > 1
        error('stats:kmeans:AmbiguousStart', ...
              'Ambiguous ''start'' parameter value:  %s.', start);
    elseif isempty(i)
        error('stats:kmeans:UnknownStart', ...
              'Unknown ''start'' parameter value:  %s.', start);
    elseif isempty(k)
        error('stats:kmeans:MissingK', ...
              'You must specify the number of clusters, K.');
    end
    start = startNames{i};
    if strcmp(start, 'uniform')
        if strcmp(distance, 'hamming')
            error('stats:kmeans:UniformStartForHamm', ...
                  'Hamming distance cannot be initialized with uniform random values.');
        end
        Xmins = min(X,[],1);
        Xmaxs = max(X,[],1);
    end
elseif isnumeric(start)
    CC = start;
    start = 'numeric';
    if isempty(k)
        k = size(CC,1);
    elseif k ~= size(CC,1);
        error('stats:kmeans:MisshapedStart', ...
              'The ''start'' matrix must have K rows.');
    elseif size(CC,2) ~= p
        error('stats:kmeans:MisshapedStart', ...
              'The ''start'' matrix must have the same number of columns as X.');
    end
    if isempty(reps)
        reps = size(CC,3);
    elseif reps ~= size(CC,3);
        error('stats:kmeans:MisshapedStart', ...
              'The third dimension of the ''start'' array must match the ''replicates'' parameter value.');
    end
    
    % Need to center explicit starting points for 'correlation'. (Re)normalization
    % for 'cosine'/'correlation' is done at each iteration.
    if isequal(distance, 'correlation')
        CC = CC - repmat(mean(CC,2),[1,p,1]);
    end
else
    error('stats:kmeans:InvalidStart', ...
          'The ''start'' parameter value must be a string or a numeric matrix or array.');
end

if ischar(emptyact)
    emptyactNames = {'error','drop','singleton'};
    i = strmatch(lower(emptyact), emptyactNames);
    if length(i) > 1
        error('stats:kmeans:AmbiguousEmptyAction', ...
              'Ambiguous ''emptyaction'' parameter value:  %s.', emptyact);
    elseif isempty(i)
        error('stats:kmeans:UnknownEmptyAction', ...
              'Unknown ''emptyaction'' parameter value:  %s.', emptyact);
    end
    emptyact = emptyactNames{i};
else
    error('stats:kmeans:InvalidEmptyAction', ...
          'The ''emptyaction'' parameter value must be a string.');
end

if ischar(display)
    i = strmatch(lower(display), strvcat('off','notify','final','iter'));
    if length(i) > 1
        error('stats:kmeans:AmbiguousDisplay', ...
              'Ambiguous ''display'' parameter value:  %s.', display);
    elseif isempty(i)
        error('stats:kmeans:UnknownDisplay', ...
              'Unknown ''display'' parameter value:  %s.', display);
    end
    display = i-1;
else
    error('stats:kmeans:InvalidDisplay', ...
          'The ''display'' parameter value must be a string.');
end

if k == 1
    error('stats:kmeans:OneCluster', ...
          'The number of clusters must be greater than 1.');
elseif n < k
    error('stats:kmeans:TooManyClusters', ...
          'X must have more rows than the number of clusters.');
end

% Assume one replicate
if isempty(reps)
    reps = 1;
end

%
% Done with input argument processing, begin clustering
%

dispfmt = '%6d\t%6d\t%8d\t%12g';
D = repmat(NaN,n,k);   % point-to-cluster distances
Del = repmat(NaN,n,k); % reassignment criterion
m = zeros(k,1);

totsumDBest = Inf;
for rep = 1:reps
    switch start
    case 'uniform'
        C = unifrnd(Xmins(ones(k,1),:), Xmaxs(ones(k,1),:));
        % For 'cosine' and 'correlation', these are uniform inside a subset
        % of the unit hypersphere.  Still need to center them for
        % 'correlation'.  (Re)normalization for 'cosine'/'correlation' is
        % done at each iteration.
        if isequal(distance, 'correlation')
            C = C - repmat(mean(C,2),1,p);
        end
        if isa(X,'single')
            C = single(C);
        end
    case 'sample'
        C = X(randsample(n,k),:);
        if ~isfloat(C)      % X may be logical
            C = double(C);
        end
    case 'cluster'
        Xsubset = X(randsample(n,floor(.1*n)),:);
        [dum, C] = mykmeans(Xsubset, k, varargin{:}, 'start','sample', 'replicates',1);
    case 'numeric'
        C = CC(:,:,rep);
    end    
    changed = 1:k; % everything is newly assigned
    idx = zeros(n,1);
    totsumD = Inf;
    
    if display > 2 % 'iter'
        disp(sprintf('  iter\t phase\t     num\t         sum'));
    end
    
    %
    % Begin phase one:  batch reassignments
    %
    
    converged = false;
    iter = 0;
    while true
        % Compute the distance from every point to each cluster centroid
        D(:,changed) = distfun(X, C(changed,:), distance, iter);
        
        % Compute the total sum of distances for the current configuration.
        % Can't do it first time through, there's no configuration yet.
        if iter > 0
            totsumD = sum(D((idx-1)*n + (1:n)'));
            % Test for a cycle: if objective is not decreased, back out
            % the last step and move on to the single update phase
            if prevtotsumD <= totsumD
                idx = previdx;
                [C(changed,:), m(changed)] = gcentroids(X, idx, changed, distance, Xsort, Xord);
                iter = iter - 1;
                break;
            end
            if display > 2 % 'iter'
                disp(sprintf(dispfmt,iter,1,length(moved),totsumD));
            end
            if iter >= maxit, break; end
        end

        % Determine closest cluster for each point and reassign points to clusters
        previdx = idx;
        prevtotsumD = totsumD;
        [d, nidx] = min(D, [], 2);

        if iter == 0
            % Every point moved, every cluster will need an update
            moved = 1:n;
            idx = nidx;
            changed = 1:k;
        else
            % Determine which points moved
            moved = find(nidx ~= previdx);
            if length(moved) > 0
                % Resolve ties in favor of not moving
                moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
            end
            if length(moved) == 0
                break;
            end
            idx(moved) = nidx(moved);

            % Find clusters that gained or lost members
            changed = unique([idx(moved); previdx(moved)])';
        end

        % Calculate the new cluster centroids and counts.
        [C(changed,:), m(changed)] = gcentroids(X, idx, changed, distance, Xsort, Xord);
        iter = iter + 1;
        
        % Deal with clusters that have just lost all their members
        empties = changed(m(changed) == 0);
        if ~isempty(empties)
            switch emptyact
            case 'error'
                error('stats:kmeans:EmptyCluster', ...
                      'Empty cluster created at iteration %d.',iter);
            case 'drop'
                % Remove the empty cluster from any further processing
                D(:,empties) = NaN;
                changed = changed(m(changed) > 0);
                if display > 0
                    warning('stats:kmeans:EmptyCluster', ...
                            'Empty cluster created at iteration %d.',iter);
                end
            case 'singleton'
                if display > 0
                    warning('stats:kmeans:EmptyCluster', ...
                            'Empty cluster created at iteration %d.',iter);
                end
                
                for i = empties
                    % Find the point furthest away from its current cluster.
                    % Take that point out of its cluster and use it to create
                    % a new singleton cluster to replace the empty one.
                    [dlarge, lonely] = max(d);
                    from = idx(lonely); % taking from this cluster
                    C(i,:) = X(lonely,:);
                    m(i) = 1;
                    idx(lonely) = i;
                    d(lonely) = 0;
                    
                    % Update clusters from which points are taken
                    [C(from,:), m(from)] = gcentroids(X, idx, from, distance, Xsort, Xord);
                    changed = unique([changed from]);
                end
            end
        end
    end % phase one

    % Initialize some cluster information prior to phase two
    switch distance
    case 'cityblock'
        Xmid = zeros([k,p,2]);
        for i = 1:k
            if m(i) > 0
                % Separate out sorted coords for points in i'th cluster,
                % and save values above and below median, component-wise
                Xsorted = reshape(Xsort(idx(Xord)==i), m(i), p);
                nn = floor(.5*m(i));
                if mod(m(i),2) == 0
                    Xmid(i,:,1:2) = Xsorted([nn, nn+1],:)';
                elseif m(i) > 1
                    Xmid(i,:,1:2) = Xsorted([nn, nn+2],:)';
                else
                    Xmid(i,:,1:2) = Xsorted([1, 1],:)';
                end
            end
        end
    case 'hamming'
        Xsum = zeros(k,p);
        for i = 1:k
            if m(i) > 0
                % Sum coords for points in i'th cluster, component-wise
                Xsum(i,:) = sum(X(idx==i,:), 1);
            end
        end
    end
    
    %
    % Begin phase two:  single reassignments
    %
    changed = find(m' > 0);
    lastmoved = 0;
    nummoved = 0;
    iter1 = iter;
    while iter < maxit
        % Calculate distances to each cluster from each point, and the
        % potential change in total sum of errors for adding or removing
        % each point from each cluster.  Clusters that have not changed
        % membership need not be updated.
        %
        % Singleton clusters are a special case for the sum of dists
        % calculation.  Removing their only point is never best, so the
        % reassignment criterion had better guarantee that a singleton
        % point will stay in its own cluster.  Happily, we get
        % Del(i,idx(i)) == 0 automatically for them.
		switch distance
		case 'sqeuclidean'
            for i = changed
                mbrs = (idx == i);
                sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
                if m(i) == 1
                    sgn(mbrs) = 0; % prevent divide-by-zero for singleton mbrs
                end
                Del(:,i) = (m(i) ./ (m(i) + sgn)) .* sum((X - C(repmat(i,n,1),:)).^2, 2);
            end
        case 'cityblock'
            for i = changed
                if mod(m(i),2) == 0 % this will never catch singleton clusters
                    ldist = Xmid(repmat(i,n,1),:,1) - X;
                    rdist = X - Xmid(repmat(i,n,1),:,2);
                    mbrs = (idx == i);
                    sgn = repmat(1-2*mbrs, 1, p); % -1 for members, 1 for nonmembers
                    Del(:,i) = sum(max(0, max(sgn.*rdist, sgn.*ldist)), 2);
                else
                    Del(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
                end
            end
        case {'cosine','correlation'}
            % The points are normalized, centroids are not, so normalize them
            normC(changed) = sqrt(sum(C(changed,:).^2, 2));
            if any(normC < eps(class(normC))) % small relative to unit-length data points
                error('stats:kmeans:ZeroCentroid', ...
                      'Zero cluster centroid created at iteration %d.',iter);
            end
            % This can be done without a loop, but the loop saves memory allocations
            for i = changed
                XCi = X * C(i,:)';
                mbrs = (idx == i);
                sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
                Del(:,i) = 1 + sgn .*...
                      (m(i).*normC(i) - sqrt((m(i).*normC(i)).^2 + 2.*sgn.*m(i).*XCi + 1));
            end
        case 'hamming'
            for i = changed
                if mod(m(i),2) == 0 % this will never catch singleton clusters
                    % coords with an unequal number of 0s and 1s have a
                    % different contribution than coords with an equal
                    % number
                    unequal01 = find(2*Xsum(i,:) ~= m(i));
                    numequal01 = p - length(unequal01);
                    mbrs = (idx == i);
                    Di = abs(X(:,unequal01) - C(repmat(i,n,1),unequal01));
                    Del(:,i) = (sum(Di, 2) + mbrs*numequal01) / p;
                else
                    Del(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
                end
            end
		end

        % Determine best possible move, if any, for each point.  Next we
        % will pick one from those that actually did move.
        previdx = idx;
        prevtotsumD = totsumD;
        [minDel, nidx] = min(Del, [], 2);
        moved = find(previdx ~= nidx);
        if length(moved) > 0
            % Resolve ties in favor of not moving
            moved = moved(Del((previdx(moved)-1)*n + moved) > minDel(moved));
        end
        if length(moved) == 0
            % Count an iteration if phase 2 did nothing at all, or if we're
            % in the middle of a pass through all the points
            if (iter - iter1) == 0 | nummoved > 0
                iter = iter + 1;
                if display > 2 % 'iter'
                    disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
                end
            end
            converged = true;
            break;
        end
        
        % Pick the next move in cyclic order
        moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;
        
        % If we've gone once through all the points, that's an iteration
        if moved <= lastmoved
            iter = iter + 1;
            if display > 2 % 'iter'
                disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
            end
            if iter >= maxit, break; end
            nummoved = 0;
        end
        nummoved = nummoved + 1;
        lastmoved = moved;
        
        oidx = idx(moved);
        nidx = nidx(moved);
        totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);
        
        % Update the cluster index vector, and rhe old and new cluster
        % counts and centroids
        idx(moved) = nidx;
        m(nidx) = m(nidx) + 1;
        m(oidx) = m(oidx) - 1;
        switch distance
        case 'sqeuclidean'
            C(nidx,:) = C(nidx,:) + (X(moved,:) - C(nidx,:)) / m(nidx);
            C(oidx,:) = C(oidx,:) - (X(moved,:) - C(oidx,:)) / m(oidx);
        case 'cityblock'
            for i = [oidx nidx]
                % Separate out sorted coords for points in each cluster.
                % New centroid is the coord median, save values above and
                % below median.  All done component-wise.
                Xsorted = reshape(Xsort(idx(Xord)==i), m(i), p);
                nn = floor(.5*m(i));
                if mod(m(i),2) == 0
                    C(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
                    Xmid(i,:,1:2) = Xsorted([nn, nn+1],:)';
                else
                    C(i,:) = Xsorted(nn+1,:);
                    if m(i) > 1
                        Xmid(i,:,1:2) = Xsorted([nn, nn+2],:)';
                    else
                        Xmid(i,:,1:2) = Xsorted([1, 1],:)';
                    end
                end
            end
        case {'cosine','correlation'}
            C(nidx,:) = C(nidx,:) + (X(moved,:) - C(nidx,:)) / m(nidx);
            C(oidx,:) = C(oidx,:) - (X(moved,:) - C(oidx,:)) / m(oidx);
        case 'hamming'
            % Update summed coords for points in each cluster.  New
            % centroid is the coord median.  All done component-wise.
            Xsum(nidx,:) = Xsum(nidx,:) + X(moved,:);
            Xsum(oidx,:) = Xsum(oidx,:) - X(moved,:);
            C(nidx,:) = .5*sign(2*Xsum(nidx,:) - m(nidx)) + .5;
            C(oidx,:) = .5*sign(2*Xsum(oidx,:) - m(oidx)) + .5;
        end
        changed = sort([oidx nidx]);
    end % phase two
    
    if (~converged) & (display > 0)
        warning('stats:kmeans:FailedToConverge', ...
                'Failed to converge in %d iterations.', maxit);
    end

    % Calculate cluster-wise sums of distances
    nonempties = find(m(:)'>0);
    D(:,nonempties) = distfun(X, C(nonempties,:), distance, iter);
    d = D((idx-1)*n + (1:n)');
    sumD = zeros(k,1);
    for i = 1:k
        sumD(i) = sum(d(idx == i));
    end
    if display > 1 % 'final' or 'iter'
        disp(sprintf('%d iterations, total sum of distances = %g',iter,totsumD));
    end

    % Save the best solution so far
    if totsumD < totsumDBest
        totsumDBest = totsumD;
        idxBest = idx;
        Cbest = C;
        sumDBest = sumD;
        if nargout > 3
            Dbest = D;
        end
    end
end

% Return the best solution
idx = idxBest;
C = Cbest;
sumD = sumDBest;
if nargout > 3
    D = Dbest;
end


%------------------------------------------------------------------

function D = distfun(X, C, dist, iter)
%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1));
nclusts = size(C,1);

switch dist
case 'sqeuclidean'
    for i = 1:nclusts
        D(:,i) = sum((X - C(repmat(i,n,1),:)).^2, 2);
    end
case 'cityblock'
    for i = 1:nclusts
        D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
    end
case {'cosine','correlation'}
    % The points are normalized, centroids are not, so normalize them
    normC = sqrt(sum(C.^2, 2));
    if any(normC < eps(class(normC))) % small relative to unit-length data points
        error('stats:kmeans:ZeroCentroid', ...
              'Zero cluster centroid created at iteration %d.',iter);
    end
    % This can be done without a loop, but the loop saves memory allocations
    for i = 1:nclusts
        D(:,i) = 1 - (X * C(i,:)') ./ normC(i);
    end
case 'hamming'
    for i = 1:nclusts
        D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
    end
end


%------------------------------------------------------------------

function [centroids, counts] = gcentroids(X, index, clusts, dist, Xsort, Xord)
%GCENTROIDS Centroids and counts stratified by group.
[n,p] = size(X);
num = length(clusts);
centroids = repmat(NaN, [num p]);
counts = zeros(num,1);
for i = 1:num
    members = find(index == clusts(i));
    if length(members) > 0
        counts(i) = length(members);
        switch dist
        case 'sqeuclidean'
            centroids(i,:) = sum(X(members,:),1) / counts(i);
        case 'cityblock'
            % Separate out sorted coords for points in i'th cluster,
            % and use to compute a fast median, component-wise
            Xsorted = reshape(Xsort(index(Xord)==clusts(i)), counts(i), p);
            nn = floor(.5*counts(i));
            if mod(counts(i),2) == 0
                centroids(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
            else
                centroids(i,:) = Xsorted(nn+1,:);
            end
        case {'cosine','correlation'}
            centroids(i,:) = sum(X(members,:),1) / counts(i); % unnormalized
        case 'hamming'
            % Compute a fast median for binary data, component-wise
            centroids(i,:) = .5*sign(2*sum(X(members,:), 1) - counts(i)) + .5;
        end
    end
end


%------------------------------------------------------------------

function [eid,emsg,varargout]=statgetargs(pnames,dflts,varargin) 
 %STATGETARGS Process parameter name/value pairs for statistics functions 
 %   [EID,EMSG,A,B,...]=STATGETARGS(PNAMES,DFLTS,'NAME1',VAL1,'NAME2',VAL2,...) 
 %   accepts a cell array PNAMES of valid parameter names, a cell array 
 %   DFLTS of default values for the parameters named in PNAMES, and 
 %   additional parameter name/value pairs.  Returns parameter values A,B,... 
 %   in the same order as the names in PNAMES.  Outputs corresponding to 
 %   entries in PNAMES that are not specified in the name/value pairs are 
 %   set to the corresponding value from DFLTS.  If nargout is equal to 
 %   length(PNAMES)+1, then unrecognized name/value pairs are an error.  If 
 %   nargout is equal to length(PNAMES)+2, then all unrecognized name/value 
 %   pairs are returned in a single cell array following any other outputs. 
 % 
 %   EID and EMSG are empty if the arguments are valid.  If an error occurs, 
 %   EMSG is the text of an error message and EID is the final component 
 %   of an error message id.  STATGETARGS does not actually throw any errors, 
 %   but rather returns EID and EMSG so that the caller may throw the error. 
 %   Outputs will be partially processed after an error occurs. 
 % 
 %   This utility is used by some Statistics Toolbox functions to process 
 %   name/value pair arguments. 
 % 
 %   Example: 
 %       pnames = {'color' 'linestyle', 'linewidth'} 
 %       dflts  = {    'r'         '_'          '1'} 
 %       varargin = {{'linew' 2 'nonesuch' [1 2 3] 'linestyle' ':'} 
 %       [eid,emsg,c,ls,lw] = statgetargs(pnames,dflts,varargin{:})    % error 
 %       [eid,emsg,c,ls,lw,ur] = statgetargs(pnames,dflts,varargin{:}) % ok 
  
 %   Copyright 1993-2004 The MathWorks, Inc.  
 %   $Revision: 1.4.2.1 $  $Date: 2003/11/01 04:28:41 $  
  
 % We always create (nparams+2) outputs: 
 %    one each for emsg and eid 
 %    nparams varargs for values corresponding to names in pnames 
 % If they ask for one more (nargout == nparams+3), it's for unrecognized 
 % names/values 
  
 % Initialize some variables 
 emsg = ''; 
 eid = ''; 
 nparams = length(pnames); 
 varargout = dflts; 
 unrecog = {}; 
 nargs = length(varargin); 
  
 % Must have name/value pairs 
 if mod(nargs,2)~=0 
     eid = 'WrongNumberArgs'; 
     emsg = 'Wrong number of arguments.'; 
 else 
     % Process name/value pairs 
     for j=1:2:nargs 
         pname = varargin{j}; 
         if ~ischar(pname) 
             eid = 'BadParamName'; 
             emsg = 'Parameter name must be text.'; 
             break; 
         end 
         i = strmatch(lower(pname),pnames); 
         if isempty(i) 
             % if they've asked to get back unrecognized names/values, add this 
             % one to the list 
             if nargout > nparams+2 
                 unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}}; 
                  
                 % otherwise, it's an error 
             else 
                 eid = 'BadParamName'; 
                 emsg = sprintf('Invalid parameter name:  %s.',pname); 
                 break; 
             end 
         elseif length(i)>1 
             eid = 'BadParamName'; 
             emsg = sprintf('Ambiguous parameter name:  %s.',pname); 
             break; 
         else 
             varargout{i} = varargin{j+1}; 
         end 
     end 
 end 
  
 varargout{nparams+1} = unrecog; 


 
 
%------------------------------------------------------------------
 
 function y = randsample(n, k, replace, w)
%RANDSAMPLE Random sample, with or without replacement.
%   Y = RANDSAMPLE(N,K) returns Y as a vector of values sampled uniformly
%   at random, without replacement, from the integers 1:N.
%
%   Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at
%   random, without replacement, from the values in the vector POPULATION.
%
%   Y = RANDSAMPLE(...,REPLACE) returns a sample taken with replacement if
%   REPLACE is true, or without replacement if REPLACE is false (the default).
%
%   Y = RANDSAMPLE(...,true,W) returns a weighted sample, using positive
%   weights W, taken with replacement.  W is often a vector of probabilities.
%   This function does not support weighted sampling without replacement.
%
%   Example:  Generate a random sequence of the characters ACGT, with
%   replacement, according to specified probabilities.
%
%      R = randsample('ACGT',48,true,[0.15 0.35 0.35 0.15])
%
%   See also RAND, RANDPERM.

%   Copyright 1993-2007 The MathWorks, Inc.
%   $Revision: 1.1.4.2 $  $Date: 2007/08/03 21:42:47 $

if nargin < 2
    error('stats:randsample:TooFewInputs','Requires two input arguments.');
elseif numel(n) == 1
    population = [];
else
    population = n;
    n = numel(population);
    if length(population)~=n
       error('stats:randsample:BadPopulation','POPULATION must be a vector.');
    end
end

if nargin < 3
    replace = false;
end

if nargin < 4
    w = [];
elseif ~isempty(w)
    if length(w) ~= n
        if isempty(population)
            error('stats:randsample:InputSizeMismatch',...
                  'W must have length equal to N.');
        else
            error('stats:randsample:InputSizeMismatch',...
                  'W must have the same length as the population.');
        end
    else
        p = w(:)' / sum(w);
    end
end

switch replace
   
% Sample with replacement
case {true, 'true', 1}
    if isempty(w)
        y = ceil(n .* rand(k,1));
    else
        [dum, y] = histc(rand(k,1),[0 cumsum(p)]);
    end

% Sample without replacement
case {false, 'false', 0}
    if k > n
        if isempty(population)
            error('stats:randsample:SampleTooLarge',...
        'K must be less than or equal to N for sampling without replacement.');
        else
            error('stats:randsample:SampleTooLarge',...
                  'K must be less than or equal to the population size.');
        end
    end
    
    if isempty(w)
        % If the sample is a sizeable fraction of the population,
        % just randomize the whole population (which involves a full
        % sort of n random values), and take the first k.
        if 4*k > n
            rp = randperm(n);
            y = rp(1:k);
            
        % If the sample is a small fraction of the population, a full sort
        % is wasteful.  Repeatedly sample with replacement until there are
        % k unique values.
        else
            x = zeros(1,n); % flags
            sumx = 0;
            while sumx < k
                x(ceil(n * rand(1,k-sumx))) = 1; % sample w/replacement
                sumx = sum(x); % count how many unique elements so far
            end
            y = find(x > 0);
            y = y(randperm(k));
        end
    else
        error('stats:randsample:NoWeighting',...
              'Weighted sampling without replacement is not supported.');
    end
otherwise
    error('stats:randsample:BadReplaceValue',...
          'REPLACE must be either true or false.');
end

if ~isempty(population)
    y = population(y);
else
    y = y(:);
end

function y = mvnpdf(X, Mu, Sigma)
%MVNPDF Multivariate normal probability density function (pdf).
%   Y = MVNPDF(X) returns the probability density of the multivariate normal
%   distribution with zero mean and identity covariance matrix, evaluated at
%   each row of X.  Rows of the N-by-D matrix X correspond to observations or
%   points, and columns correspond to variables or coordinates.  Y is an
%   N-by-1 vector.
%
%   Y = MVNPDF(X,MU) returns the density of the multivariate normal
%   distribution with mean MU and identity covariance matrix, evaluated
%   at each row of X.  MU is a 1-by-D vector, or an N-by-D matrix, in which
%   case the density is evaluated for each row of X with the corresponding
%   row of MU.  MU can also be a scalar value, which MVNPDF replicates to
%   match the size of X.
%
%   Y = MVNPDF(X,MU,SIGMA) returns the density of the multivariate normal
%   distribution with mean MU and covariance SIGMA, evaluated at each row
%   of X.  SIGMA is a D-by-D matrix, or an D-by-D-by-N array, in which case
%   the density is evaluated for each row of X with the corresponding page
%   of SIGMA, i.e., MVNPDF computes Y(I) using X(I,:) and SIGMA(:,:,I).
%   Pass in the empty matrix for MU to use its default value when you want
%   to only specify SIGMA.
%
%   If X is a 1-by-D vector, MVNPDF replicates it to match the leading
%   dimension of MU or the trailing dimension of SIGMA.
%
%   Example:
%
%      mu = [1 -1]; Sigma = [.9 .4; .4 .3];
%      [X1,X2] = meshgrid(linspace(-1,3,25)', linspace(-3,1,25)');
%      X = [X1(:) X2(:)];
%      p = mvnpdf(X, mu, Sigma);
%      surf(X1,X2,reshape(p,25,25));
%
%   See also MVTPDF, MVNCDF, MVNRND, NORMPDF.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.2.4.7 $  $Date: 2006/11/11 22:55:31 $

if nargin<1
    error('stats:mvnpdf:TooFewInputs','Requires at least one input.');
elseif ndims(X)~=2
    error('stats:mvnpdf:InvalidData','X must be a matrix.');
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[n,d] = size(X);
if d<1
    error('stats:mvnpdf:TooFewDimensions','X must have at least one column.');
end

% Assume zero mean, data are already centered
if nargin < 2 || isempty(Mu)
    X0 = X;

% Get scalar mean, and use it to center data
elseif numel(Mu) == 1
    X0 = X - Mu;

% Get vector mean, and use it to center data
elseif ndims(Mu) == 2
    [n2,d2] = size(Mu);
    if d2 ~= d % has to have same number of coords as X
        error('stats:mvnpdf:InputSizeMismatch',...
              'X and MU must have the same number of columns.');
    elseif n2 == n % lengths match
        X0 = X - Mu;
    elseif n2 == 1 % mean is a single row, rep it out to match data
        X0 = X - repmat(Mu,n,1);
    elseif n == 1 % data is a single row, rep it out to match mean
        n = n2;
        X0 = repmat(X,n2,1) - Mu;
    else % sizes don't match
        error('stats:mvnpdf:InputSizeMismatch',...
              'X or MU must be a row vector, or X and MU must have the same number of rows.');
    end
    
else
    error('stats:mvnpdf:BadMu','MU must be a matrix.');
end

% Assume identity covariance, data are already standardized
if nargin < 3 || isempty(Sigma)
    % Special case: if Sigma isn't supplied, then interpret X
    % and Mu as row vectors if they were both column vectors
    if (d == 1) && (numel(X) > 1)
        X0 = X0';
        [n,d] = size(X0);
    end
    xRinv = X0;
    logSqrtDetSigma = 0;
    
% Single covariance matrix
elseif ndims(Sigma) == 2
    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and Mu as row vectors if they were both column vectors.
    if (d == 1) && (numel(X) > 1) && (size(Sigma,1) == n)
        X0 = X0';
        [n,d] = size(X0);
    end
    
    % Make sure Sigma is the right size
    sz = size(Sigma);
    if sz(1) ~= sz(2)
        error('stats:mvnpdf:BadCovariance',...
              'SIGMA must be a square matrix.');
    elseif ~isequal(sz, [d d])
        error('stats:mvnpdf:InputSizeMismatch',...
              'SIGMA must be a square matrix with size equal to the number of columns in X.');
    else
        % Make sure Sigma is a valid covariance matrix
        [R,err] = cholcov(Sigma,0);
        if err ~= 0
            error('stats:mvnpdf:BadCovariance', ...
                  'SIGMA must be symmetric and positive definite.');
        end
        % Create array of standardized data, and compute log(sqrt(det(Sigma)))
        xRinv = X0 / R;
        logSqrtDetSigma = sum(log(diag(R)));
    end
    
% Multiple covariance matrices
elseif ndims(Sigma) == 3
    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and Mu as row vectors if they were both column vectors.
    if (d == 1) && (numel(X) > 1) && (size(Sigma,1) == n)
        X0 = X0';
        [n,d] = size(X0);
    end
    
    % Data and mean are a single row, rep them out to match covariance
    if n == 1 % already know size(Sigma,3) > 1
        n = size(Sigma,3);
        X0 = repmat(X0,n,1); % rep centered data out to match cov
    end

    % Make sure Sigma is the right size
    sz = size(Sigma);
    if sz(1) ~= sz(2)
        error('stats:mvnpdf:BadCovariance',...
              'Each page of SIGMA must be a square matrix.');
    elseif (sz(1) ~= d) || (sz(2) ~= d) % Sigma is a stack of dxd matrices
        error('stats:mvnpdf:InputSizeMismatch',...
              'Each page of SIGMA must be a square matrix with size equal to the number of columns in X.');
    elseif sz(3) ~= n
        error('stats:mvnpdf:InputSizeMismatch',...
              'SIGMA must have one page for each row of X.');
    else
        
        % Create array of standardized data, and vector of log(sqrt(det(Sigma)))
        xRinv = zeros(n,d,superiorfloat(X0,Sigma));
        logSqrtDetSigma = zeros(n,1,class(Sigma));
        for i = 1:n
            % Make sure Sigma is a valid covariance matrix
            [R,err] = cholcov(Sigma(:,:,i),0);
            if err ~= 0
                error('stats:mvnpdf:BadCovariance',...
                      'Each page of SIGMA must be symmetric and positive definite.');
            end
            xRinv(i,:) = X0(i,:) / R;
            logSqrtDetSigma(i) = sum(log(diag(R)));
        end
    end
   
elseif ndims(Sigma) > 3
    error('stats:mvnpdf:BadCovariance',...
          'SIGMA must be a matrix or a 3 dimensional array.');
end

% The quadratic form is the inner products of the standardized data
quadform = sum(xRinv.^2, 2);

y = exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);

%-----------------------------------------------------------------------

function [T,p] = cholcov(Sigma,flag)
%CHOLCOV  Cholesky-like decomposition for covariance matrix.
%   T = CHOLCOV(SIGMA) computes T such that SIGMA = T'*T.  SIGMA must be
%   square, symmetric, and positive semi-definite.  If SIGMA is positive
%   definite, then T is the square, upper triangular Cholesky factor.
%
%   If SIGMA is not positive definite, T is computed from an eigenvalue
%   decomposition of SIGMA.  T is not necessarily triangular or square in
%   this case.  Any eigenvectors whose corresponding eigenvalue is close to
%   zero (within a small tolerance) are omitted.  If any remaining
%   eigenvectors are negative, T is empty.
%
%   [T,P] = CHOLCOV(SIGMA) returns the number of negative eigenvalues of
%   SIGMA, and T is empty if P>0.  If P==0, SIGMA is positive semi-definite.
%
%   If SIGMA is not square and symmetric, P is NaN and T is empty.
%
%   [T,P] = CHOLCOV(SIGMA,0) returns P==0 if SIGMA is positive definite, and
%   T is the Cholesky factor.  If SIGMA is not positive definite, P is a
%   positive integer and T is empty.  [...] = CHOLCOV(SIGMA,1) is equivalent
%   to [...] = CHOLCOV(SIGMA).
%
%   Example:
%   Factor a rank-deficient covariance matrix C.
%       C = [2 1 1 2;1 2 1 2;1 1 2 2;2 2 2 3]
%       T = cholcov(C)
%       C2 = T'*T
%   Generate data with this covariance (aside from random variation).
%       C3 = cov(randn(10000,3)*T)
%
%   See also CHOL.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2006/12/15 19:29:55 $

if nargin < 2, flag = 1; end

% Test for square, symmetric
[n,m] = size(Sigma);
wassparse = issparse(Sigma);
tol = 10*eps(full(max(abs(diag(Sigma)))));
if (n == m) && all(all(abs(Sigma - Sigma') < tol))
    [T,p] = chol(Sigma);

    % Test for positive definiteness
    if p > 0 && flag
        % Can get factors of the form Sigma==T'*T using the eigenvalue
        % decomposition of a symmetric matrix, so long as the matrix
        % is positive semi-definite.
        [U,D] = eig(full((Sigma+Sigma')/2));

        % Pick eigenvector direction so max abs coordinate is positive
        [ignore,maxind] = max(abs(U),[],1);
        negloc = (U(maxind + (0:n:(m-1)*n)) < 0);
        U(:,negloc) = -U(:,negloc);

        D = diag(D);
        tol = eps(max(D)) * length(D);
        t = (abs(D) > tol);
        D = D(t);
        p = sum(D<0); % number of negative eigenvalues

        if (p==0)
            T = diag(sqrt(D)) * U(:,t)';
        else
            T = [];
        end
    end

else
    T = [];
    p = NaN;
end

if wassparse
    T = sparse(T);
end
