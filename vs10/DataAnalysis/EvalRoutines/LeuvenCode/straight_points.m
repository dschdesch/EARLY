% This function re-aligns the data obtained of a corner shaped object into 
% a perfectly alligned object.
% This is necessary in order to get an average profile of the object with
% minimum error or noise and it also simplifies the laboratory work.
% Please note this function only works with small angles. The greater the 
% angle the larger the error.

function [x,yRot,zRot]= straight_points(x,y,z)


% Making data positive
min_z = min(z);
z = z + abs(min_z);   % z points are all positive now

% Store Matrix in order
M = [x y z];
M = sortrows(M);          
j=1;g=0;h=0;i=1;

% Returns size of each y profile
    while (g == h)
        g = M(i,j);
        h = M(i+1,j);
        i = i+1;
    end
      i = i-1; 
      
% Storing information by profiles
m = length(M);        
Profile = mat2cell(M,repmat(i,[m/i 1]),3);

% Check Highest Point (Global Maximum)
ascending = sortrows(M,3);
descending = flipdim(ascending,1);
highest_p = descending(1,3);

% Check Highest point of each profile(Local Maximums)
n=752;
asc{n}  = 0;         % Pre - allocating
x_mp{n} = 0;         %   variables            (FOR SPEED)
y_mp{n} = 0;         %   :D
z_mp{n} = 0;
p = length(Profile);
for n=1:p
 asc{n} = flipdim(sortrows(Profile{n},3),1);

 x_mp{n} = (asc{n}(1,1));
 y_mp{n} = (asc{n}(1,2));
 z_mp{n} = (asc{n}(1,3));

end

 max_xp = (cell2mat(x_mp))';
 max_yp = (cell2mat(y_mp))';
 max_zp = (cell2mat(z_mp))';
 
 maximum_p =[max_xp max_yp max_zp];


% Reference Line through it
% XZ Plane
line_XZ = polyfit(maximum_p(:,1) ,maximum_p(:,3) , 1);
line_XZ_data = polyval(line_XZ ,maximum_p(:,1));

slope =  line_XZ(1);         
theta = atan(slope);

% Distance from Max points to the reference line
distance_max_p = maximum_p(:,3) - line_XZ_data;

% New Reference Line
z_line_h(1:p,1) = maximum_p(1,3);

% Transforming Points
transform_max_p = distance_max_p + highest_p; 
    
% Rotating Data
zRot = x*sin(-theta) + z*cos(-theta);


% Reference Line
% XY Plane
%
line_XY = polyfit(max_xp, max_yp, 1);
line_XY_data = polyval(line_XY,x); 

slope2 =  line_XY(1);         
theta2 = atan(slope2);

% Rotating Data
yRot = x*sin(-theta2) + y*cos(-theta2);
yRotmax = max_xp*sin(-theta2) + max_yp*cos(-theta2);
fit_h = polyfit(max_xp,yRotmax,1);
fit_h_data = polyval(fit_h,max_xp);


%  
% % Luis Alberto Canizares