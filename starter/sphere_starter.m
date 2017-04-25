%% Starter code for source localization on spherical head models
%%
%    Data science
%    Homework 3
%    Author : Ferdinand Dowouo
%    Problem 2 : Source Localization
%%
clear; close all , clc;
warning('off')
dipole_L = 10;
[dipole_grid, thetas, phis] = create_grid(dipole_L);

sensor_L = 15;
sensor_grid = create_grid(sensor_L);
m = sensor_L^2;

A = leadfield_matrix(dipole_grid, sensor_grid);

x = zeros(dipole_L^2, 1);
x(5^2+5+1) = 1;
sigma = 0.0;
y = A * x + sigma * randn(m, 1);
figure;
plot_x_hat(x, thetas, phis);
title('Original spherical model for dipole = 10');
%% Code for l2-regularized least squares inverse solution

%dipole_L = ...;                    % Reset dipole_L to create the source space


dipole = [5, 10 , 15] ;

for d = 1: length(dipole)
    dipole_L = dipole(d) ;
    [dipole_grid, thetas, phis] = create_grid(dipole_L);
    sensor_L = 15;
    sensor_grid = create_grid(sensor_L);
    m = sensor_L^2;
    A = leadfield_matrix(dipole_grid, sensor_grid);
%     x = zeros(dipole_L^2, 1);
%     x(5^2+5+1) = 1;
%     sigma = 0.0;
%     y = A * x + sigma * randn(m, 1);
  % now entering the loop to calculate the error 
  lambda = logspace(-6, -2 , 20) ;
  numIter = 10 ;
  Error = zeros(length(lambda) ,1);
  for lbd = 1 : length(lambda)
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
              %             Iden = eye(length(train_data)) ;
              Iden = eye(size(Atrain,2)) ;
              
              W = ((Atrain')*Atrain + lambda(lbd)*Iden  ) \ (Atrain' *train_data) ;
              err = err+  sum((test_data - Atest*W).^2) ;
          end
          localError(iter) = err/nFold ;
      end
      Error(lbd)=mean(localError) ;
  end
  
  
  [minError , index] = min(Error) ;
  
  % the optimum lambda is therefore the lambda that correspond to the index
  opt_lambda = lambda(index) ;
  fprintf( 'The optimum lambda for this regularization is : %f \n ' , opt_lambda)
  
  figure;
  plot( (1:length(lambda)), Error , 'r', 'LineWidth' , 1.3) ;
  xlim([1,length(lambda)]) ;
  title(strcat('Error vs Lambda for dipole = ' , num2str(dipole_L), ' and lambda opt = ' , num2str(opt_lambda))) ;
  W_opt = ((Atrain')*Atrain + opt_lambda*Iden  ) \ (Atrain' *train_data) ;
  
  % Plot reconstruction assuming solution is stored in x_hat
  figure;
  plot_x_hat(W_opt, thetas, phis);
  title(strcat('Recovered spherical model for dipole = ' , num2str(dipole_L)));
end
%% Comment
% We can see that as the dipole value increase , we are getting a better
% and better recovery of the spherical head model, however we can also
% realize that the value of lambda that minimize the error is a very small
% value in the order of 1e-6. and when we plot the error vs lambda we
% notice that the error is icreasing as the value of lambda increase. 
