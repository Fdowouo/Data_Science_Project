%% Simple one-dimensional T-SVD example
%%
%    Data science
%    Homework 3
%    Author : Ferdinand Dowouo
%    Problem 1 : Understanding the SVD
%%
clear;  close all ; clc;
warning('off') ;
% Blur kernel
kernel_width = 21; % Must be odd
h = 0.5 + 0.5 * cos(linspace(-pi, pi, kernel_width));
% Sensing matrix - cyclic shifts of h form the rows of A.
n = 100;
h_cyclic = zeros(1, n);
kw_by_2 = floor(kernel_width / 2) + 1;
h_cyclic(1 : kw_by_2) = h(kw_by_2 : end);
h_cyclic(n-kw_by_2+2 : n) = h(1 : kw_by_2-1);
A = toeplitz(h_cyclic);
A = A(1:2:end, :);
m = size(A, 1);
% Random non-zero indices
nnz = 1;  % Number of non-zero indices
idx = datasample(1:n, nnz, 'Replace', false);
x = zeros(n, 1);
x(idx) = 1;
% Apply blur kernel
y_true = A * x;
% Add sensor noise
sigma = 0.1;
y = y_true + sigma * randn(m, 1);
% Plot original and sensed signals
figure(1); 
tx = 0:n-1;
ty = 0:2:n-1;
plot(tx, x, 'k-', 'Linewidth', 2);
hold on;
plot(ty, y, 'ro');

legend( 'Original' , 'sensed')
title('Original and sensed signals')

%% T-SVD code goes here
nTrunc = 50 ;
xRec = 0 ;
[U, E , V] = svd(A) ;
sig = diag(E) ;
d = y_true ;
for k = 1:nTrunc
    Ul = U(:,k) ;
    Vl = V(:, k) ;
    sigk = sig(k) ;
    xRec = xRec + ((Ul'*d)/sigk)*Vl ;
end
figure;
hold on
plot(tx, x, 'k-', 'Linewidth', 2);
plot(ty, y, 'ro');
plot(tx, xRec,  '-b', 'LineWidth' , 1.2) ;
hold off ;
legend( 'Original' , 'sensed', 'reconstructed')
title('Recovery clean signal without Threshold')
%% For small threshold
nThreshold =25 ;
xRecT = 0 ;
d = y_true ;
for k = 1:nThreshold
    Ul = U(:,k) ;
    Vl = V(:, k) ;
    sigk = sig(k) ;
    xRecT = xRecT + ((Ul'*d)/sigk)*Vl ;
end
figure ;
hold on
plot(tx, x, 'k-', 'Linewidth', 2);
plot(ty, y, 'ro');
plot(tx, xRecT,  '-b', 'LineWidth' , 1.2) ;
hold off
legend( 'Original' , 'sensed', 'reconstructed')
title( strcat('Recovery clean signal with small threshold = ' , num2str(nThreshold)))
%% TruncatedSVD for noizy signal
% let's find the 
% L is the identity matrix in this case
[U, E , V] = svd(A) ;
sig = diag(E) ;
nThreshold = 25;
xRec = 0 ;
for k = 1:nThreshold
    Ul = U(:,k) ;
    Vl = V(:, k) ;
    sigk = sig(k) ;
    xRec = xRec + ((Ul'*y)/sigk)*Vl ;
end
fPrime = xRec ;
Vk = V(:, nThreshold+1:end) ;
Vk_plus = ((Vk'*Vk))\(Vk') ;
c = -Vk_plus*fPrime ;
xRec_N = fPrime+ Vk*c  ;

figure ;
hold on
plot(tx, x, 'k-', 'Linewidth', 2);
plot(ty, y, 'ro');
plot(tx, xRec_N,  '-b', 'LineWidth' , 1.2) ;
hold off
legend( 'Original' , 'sensed', 'reconstructed')
title( strcat('Recovery noizy signal with small threshold = ' , num2str(nThreshold)))

% 
%% Group comment on i , ii, iii, iv, 

% i) If A is a 50 x 100 matrix and given that A is full rank, this mean that
% the rank of A is 50  , therefore the maximum number of singular values
% you can invert in absence of noise is 50 

%ii) for the Truncated SVD, using the 50th non zero singular value we can
%see that the recovered signal capture the original but still presents a
%higher variance around the true signal , this is due to the fact that some
%of the singular value are very small and when we divide a term such as
%(U*d) then the whole thing become very large and that affect the noiziness
%of the recovered signal.

% we can never get the exact original signal  because no matter how smalll
% the threshold is there will be a small value 
% iii) whe we use smaller threshold, the curve of the recovered signal
% become smooth , however we notice a decrease of the amplitude of peaks
% vallue. again this smoothness is due to the fact that smaller signular
% value have beeng removed and the contribution of the correcponding
% element of V have been turned to zero , this reduce the variability of
% the signal. at even very small value of the threshold, the curve is
% pretty much flat with a small maximum where the peak should be.

% iv)  when the noise is added, the smallest threshold ( <=15) produce the better
% result. Howerver, with the treshold more that 30, the recovered signal
% does not look anything closer to the original signal. this is in part due
% to the extra term that are getting explosive value due to very small
% value of sigma. 

%%  k fold cross validation

% using 5 fold cross-validation with 10 repetitions

nFold = 5 ;
dataParti= cvpartition(y, 'k' , nFold) ; 
lambda = 10:1:20; 
Error = [ ] ;
iteration = 10  ; % this is the number of time to run the partition
nnn = 100;
nnz = 1;  % Number of non-zero indices
idx = datasample(1:nnn, nnz, 'Replace', false);
x = zeros(nnn, 1);
x(idx) = 1;
mm = size(A, 1);

% Apply blur kernel
y_true = A * x;
% Add sensor noise
sigma = 0.1;
y = y_true + sigma * randn(mm, 1);

lambda = 10:1:20; 
Error = zeros(length(lambda),1) ;
numIter = 10 ;
for l = 1: length(lambda)
      localError =  zeros(numIter,1) ;
      for iter = 1: numIter
          nFold = 5 ;
          dataParti= cvpartition(y, 'k' , nFold) ;
          err = 0 ;
          for n = 1:nFold
              indicator = dataParti.training(n)  ;
              trainInd = find(indicator ==0) ; % index to be remove in d and in A
              train_data= y ;
              train_data(trainInd) = [] ;
              test_data = y(trainInd) ;
              Atrain = A ;
              Atrain(trainInd,:) = [] ;
              Atest = A(trainInd,:) ;
              [Utrain, Etrain, Vtrain ] = svd(Atrain) ;
              
              sigTrain = diag(Etrain) ;
              nThreshold = lambda(l);
              xTrainRec = 0 ;
              xRec = 0 ;
              for k = 1:nThreshold
                  Ul = Utrain(:,k) ;
                  Vl = Vtrain(:, k) ;
                  sigk = sigTrain(k) ;
                  xRec = xRec + ((Ul'*train_data)/sigk)*Vl ;
              end
            fPrime = xRec ;
            Vk = Vtrain(:, nThreshold+1:end) ;
            Vk_plus = ((Vk'*Vk))\(Vk') ;
            c = -Vk_plus*fPrime ;
            f = fPrime+ Vk*c  ;
           err = err + norm(test_data - (Atest*f)) ;
        end
       localError(iter) = err/nFold ;   
    end
    Error(l)=mean(localError) ;
end
[minErr , indErr] = min(Error) ;
lambda_opt = lambda(indErr) ;
figure;
plot( lambda, Error, 'r' , 'LineWidth' , 1.4)
title(strcat('regularization : Error vs Lambda , Opt_lambda = ' , num2str(lambda_opt)))
fprintf('The optimun lambda is : %d \n ' , lambda_opt)
%% Comment 
% we can see that the curve of the error is not as smooth as we wanted this
% is due to the very small size of the test set which gives rise to a
% higher variation 







