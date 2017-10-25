function ave_profile = average_profile(data_straight)
% This function finds the average profile of a 3 dimensional edge corner. Hence
% surface roughness is filtered and a smoother result is obtained.
 
% Rearranging data
M = sortrows(data_straight);   % Sorts data by profiles

% Setting counters to initial values
i=1;j=1;g=0;h=0;n=1;

% Returns size of a single profile
 while (g == h)
        g = M(i,j);
        h = M(i+1,j);
        i = i+1;
 end
        i(n) = i-1; 
        
% Separating data into all profiles        
m = length(M);        
Profile = mat2cell(M,repmat(i,[m/i 1]),3); 

% Average profile

number_of_profiles = length(Profile);

% average Z axis
sum_z = 0;sum_x =0;

average_x{i}= 0;
average_z{i}= 0; 
for r = 1:i;
    for m1 = 1:number_of_profiles;
        sum_x = sum_x + Profile{m1}(r,2);
        sum_z = sum_z + Profile{m1}(r,3);
    end
    
   average_x{r}= sum_x/number_of_profiles;
   average_z{r}= sum_z/number_of_profiles;

   sum_z = 0; 
   sum_x = 0;
end  
 
average_z = average_z';
average_x = average_x';

matrix_z = cell2mat(average_z);
matrix_x = cell2mat(average_x);

ave_profile = [matrix_x matrix_z];

% Luis Alberto Canizares