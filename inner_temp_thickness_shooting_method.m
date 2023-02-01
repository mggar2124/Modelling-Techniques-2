function [err,A] = inner_temp_thickness_shooting_method(initial_guess_thickness_1,initial_guess_thickness_2,target_temp,method)

%Function uses a shooting method to choose the optimum tile thickness
%through calculation with differing target temperatures and differing
%methods. Function uses two initial guesses to allow for further
%calculations. The specified target_temp is the maximum temperature at
%which the inner structure breaks/becomes damaged.

%an example [err,A] = inner_temp_thickness_shooting_method(0.01,0.03,120,'backward-neumann')   
 
n = 1;
%computation for the initial first guess
[~, ~, u] = Shuttle_Final(4000, 501, initial_guess_thickness_1, 41, method, false);
current_error = target_temp - max(u(:,1));
index(n)=n;
A(n) = initial_guess_thickness_1;
err(n) = current_error;


%computation for the initial second guess
n=n+1;
[~, ~, u] = Shuttle_Final(4000, 501, initial_guess_thickness_2, 41, method, false);
current_error = target_temp - max(u(:,1));
index(n)=n;
A(n) = initial_guess_thickness_2;
err(n) = current_error;

while abs(err(end)) > 0.001 % computes to the accuracy of milimeters - just a tolerance example, in real life such accuracy is unfeasible
    new_alpha = A(n) - err(n)*((A(n)-A(n-1))/(err(n)-err(n-1))); %function derived from the shooting method
    %uses previous two errors and associated thicknesses to estimate
    %the next guess for thickness
    [~, ~, u] = Shuttle_Final(4000, 501, new_alpha, 41, method, false);
    Res = target_temp - max(u(:,1)); %new error
    n = n+1;
    err(n) = Res;%adds latest error to the error array
    A(n) = new_alpha;% adds latest thickness to height array
    index(n)=n;
end

figure(1)
plot(index,err)
legend('Shooting Methods Error Iterations','Location','Northwest')
grid on
grid minor
xlabel('Iteration Number') 
ylabel('Error: Temperature in degrees C') 

figure(2)
plot(index,A)
legend('Thickness Value','Location','Northwest')
grid on
grid minor
xlabel('Iteration Number') 
ylabel('Thickness in Meters') 


