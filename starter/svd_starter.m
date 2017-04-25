%% Simple one-dimensional T-SVD example
clear;  close all ; clc;

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
hold off 
legend( 'Original' , 'sensed' ,'reconstructed')
title('Original and sensed signals')

%% T-SVD code goes here

[ U ,E , V_t] = svd(A) ;
V = V_t ;
figure ;
subplot(211)
hold on 
for i = 1:4
    plot(U(:, i) ) ;
end
hold off
title(' first 4 column of  U')
subplot(212) , hold on
for j = 1:4
    plot(V(: , j)) ;
    xlim([0,50]) ;
end
hold off
title(' first 4 column of V')
%% comment : what do you observe?
% the first 4 column of U are just  a straigth line and sine and cosine
% function which are representative of the cyclic function and also the
% typology of the matrix A , 

% phase shifted by 90 degree so the dot product is one

%%
kernelWidth = [ 11 , 21, 31, 41] ; 
figure ;
for ker = 1 : length(kernelWidth) 
     h = 0.5 + 0.5 * cos(linspace(-pi, pi, kernelWidth(ker)));
     n = 100;
     h_cyclic = zeros(1, n);
     kw_by_2 = floor(kernelWidth(ker) / 2) + 1;
     h_cyclic(1 : kw_by_2) = h(kw_by_2 : end);
     h_cyclic(n-kw_by_2+2 : n) = h(1 : kw_by_2-1);
     A = toeplitz(h_cyclic);
     A = A(1:2:end, :);
    [U, E , Vt] = svd(A) ;
    singValue = diag(E) ;
    subplot(2,2,ker) 
    plot((1:size(E,1)) , singValue , 'r', 'LineWidth' , 1.3) ; 
    title(strcat('Singular Value  for k =  ' , num2str( kernelWidth(ker)) ))
end
%% Comment 
% we can see that the decay is much faster as the value of kernel_width
% increase, what this tell us is that as the singular value tend to zero,
% the corresponding colunms of A does not contribute anything to the
% recovery of the data.  The more rapid the decrease of the singular value,
% the  less we can reliably reconstruct for a given noise level.


%%  Part C) : effect of singular value and noise
kernel_width = 21 ;
nnz = 1 ;

%% i)
% is A is a 50 x 100 matrix and given that A is full rank, this mean that
% the rank of A is 50  , therefore the maximum number of singular values
% you can invert in absence of noise is 50 

%% ii)  code the truncated SVD algorithm 

% x_hat = (sum(1/sigma_k)*V_k*U_k ' ) * y  ; for k = 1 , 2 ... r

kernel_width = 21 ;
nnz = 1 ;

h = 0.5 + 0.5 * cos(linspace(-pi, pi, kernel_width));
n = 100;
h_cyclic = zeros(1, n);
kw_by_2 = floor(kernel_width / 2) + 1;
h_cyclic(1 : kw_by_2) = h(kw_by_2 : end);
h_cyclic(n-kw_by_2+2 : n) = h(1 : kw_by_2-1);
A = toeplitz(h_cyclic);
A = A(1:2:end, :);
[U, E , Vt] = svd(A) ;
V = Vt' ;
singValue = diag(E) ;

nTrunc = 50 ;
x_hat = 0 ;
x_hat_true = 0 ;

E_hat = E(:, 1:50) ;
U_hat = U ;
V_hat = V(: , 1:50) ;
sig = diag( E_hat) ;
 xhat_true = (V_hat*inv(E_hat)*U_hat)*y_true ;
 xhat = (V_hat*inv(E_hat)*U_hat)*y  ;
% figure;
hold on ;
% plot(tx,xhat , '-k') ;
plot(tx, xhat_true, '-g') ;   hold off ;
% legend('noizy' , 'clean')
title(' truncated SVD without threshold') ;
%%  Seting a small treshold
threshold = 5 ;
diag_sig = diag(E) ;
remove_ind = find(diag_sig <=threshold) ;
r =  50 - length(remove_ind) ;
E_trunc = E(1:r , 1:r) ;
V_trunc = V(:, 1:r) ;
U_trunc = U(1:r , :) ; 
x_rec_true = (V_trunc*inv(E_trunc)*U_trunc)*y_true ;
x_rec  = (V_trunc*inv(E_trunc)*U_trunc)*y ;

figure; hold on , plot(tx,x_rec , '-k') ; plot(tx, x_rec_true, '-r') ;   hold off ;
legend('noizy' , 'clean')
title(' truncated SVD with small threshold') ;

 %% using a smaller threshold
 % we can see that using a smaller threshold does not really change the
 % result , using only the 7th highest singular value does not change the
 % result that much. however we still do not have full reconstruction and
 % this is due to the typology of the matrix A .
 


%%
num = 50 ;
xhat = 0 ;
xhat_t  = 0 ;
d = y ;
dt = y_true ;

[U, E , Vt] = svd(A) ;
V = Vt' ;
U_k = U(: , 1:num) ;
V_k = V(1:num, :) ;
E_k =  E(1:num , 1:num) ;
diag_k = diag(E_k) ;
for k = 1 : num
%     xhat = xhat + ( ( U_k(:,k)' *d)/diag_k(k))*V_k(:, k) ;
    xhat_t = xhat_t + ( ( U_k(:,k)' *dt)/diag_k(k))*V_k(:, k) ;
end

figure , 
hold on
plot(tx, xhat,  '-b') ;
plot(tx, xhat_t , '-r') ;
hold off 























