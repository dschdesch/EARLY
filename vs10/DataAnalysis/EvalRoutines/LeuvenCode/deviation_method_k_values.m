% This method is based on the separation of the function to the line
% tangent to the hills of the mountain. The code calculates the distance
% between the profile and the straight line generated. Once the function
% changes from one side of the line to the other (due to surface roughness)
% the counter will go back to zero. If this counter gets to a certain
% threshold, ie. the graph is no longer a straight line, the program will
% understand that this is the point where the values k1 and k2 start.


function deviation_method_k1_and_k2 = deviation_method_k_values(ave_profile,length_profile, sliders)

% Find Max point
global left_waterfall_deviation;
global right_waterfall_deviation;

global top_left_slider
global top_right_slider
global bottom_left_slider
global bottom_right_slider

top_left_slider = sliders(1,1);
top_right_slider = sliders(1,2);
bottom_left_slider = sliders(2,1);
bottom_right_slider = sliders(2,2);

descending = flipdim(sortrows(ave_profile,2),1);   
maxPoint = descending(1,1:2);
i=1;
    % Find max point in ave_profile
for count = 1:length_profile 
   if ave_profile(count,1:2)== maxPoint(1,1:2)
    break
   end
   i = i+1;
end

% High Reference Points 

left_high_clearance = top_left_slider;
left_high = ave_profile(i-left_high_clearance,1:2);

right_high_clearance = top_right_slider;
right_high = ave_profile(i+right_high_clearance,1:2);

        % 40 is considered to be a reasonable
        % amount of points for the separation
        % of the corner. i = 40 is x = ~100microns
        % Check graph and change if needed
        
% Low Reference Points
left_low_clearance = bottom_left_slider;
right_low_clearance = bottom_right_slider;
j=1;
for count = 1:length_profile
    if ave_profile(count,2)> ave_profile(1,2);
     break
    end
    j = j+1;
end 
left_low = ave_profile(j+left_low_clearance,1:2);

j2=1;
for count = length_profile:-1:1
    if ave_profile(count,2)> ave_profile(end,2);
     break
    end
    j2 = j2+1;
end
right_low = ave_profile(length_profile - j2 - right_low_clearance,1:2);

% Getting Points in Between Reference Points 
left_arm = ave_profile(j+left_low_clearance:i-left_high_clearance,1:2);
right_arm = ave_profile(i+right_high_clearance:length_profile-j2-right_low_clearance,1:2);


% Fitting lines through these points 
left_fit = polyfit(left_arm(:,1),left_arm(:,2),1);
left_points = ave_profile((j+left_low_clearance):(i+right_high_clearance),1:2);
left_line = polyval(left_fit,left_points(:,1));

right_fit = polyfit(right_arm(:,1),right_arm(:,2),1);
right_points = ave_profile((i-left_high_clearance):(length_profile - j2 - right_low_clearance),1:2);
right_line = polyval(right_fit,right_points(:,1));


% Find Waterfall Points
% LEFT 
% Testing for convergence
length_left_points = length(left_points);
distance_left(length_left_points) = 0;
j = 0;
count2 = length_left_points;

global devitolerance


deviation_clearance = devitolerance;


for count = length_left_points:-1:1
   count2 = count2 - 1;
    distance_left(count) = left_line(count) - left_points(count,2);       
    
    j = j+1;
    if distance_left(count) <= deviation_clearance
        break
    end 
end
left_waterfall_deviation = left_points(count2,1:2);

% RIGHT
% Testing for convergence
length_right_points = length(right_points);
distance_right(length_right_points) = 0;
count3 = 0;


for count = 1:length_right_points
    distance_right(count) = right_line(count) - right_points(count,2);
    count3 = count3+1;
    if distance_right(count) < deviation_clearance
       break 
    end
       
end
right_waterfall_deviation = right_points(count3,1:2);


%
% Intersection of tangent lines
y1 = right_line(1,1);      x1 = right_points(1,1);
y2 = right_line(end,1);    x2 = right_points(end,1);
y3 = left_line(1,1);       x3 = left_points(1,1); 
y4 = left_line(end,1);     x4 = left_points(end,1); 

xcomptop = (((x1*y2)-(y1*x2))*(x3-x4))-((x1-x2)*(x3*y4 - y3*x4));   % x vector component top
xcompbottom = ((x1-x2)*(y3 - y4))-((y1-y2)*(x3-x4));                % x vector component bottom  

ycomptop = (((x1*y2)-(y1*x2))*(y3-y4))-((y1-y2)*(x3*y4 - y3*x4));   % y vector component top
ycompbottom = ((x1-x2)*(y3 - y4))-((y1-y2)*(x3-x4));                % y vector component bottom

intersection = [(xcomptop/xcompbottom),(ycomptop/ycompbottom)];

k1 = sqrt((intersection(1,1)- left_waterfall_deviation(1,1))^2+(intersection(1,2)-left_waterfall_deviation(1,2))^2);
k2 = sqrt((intersection(1,1)- right_waterfall_deviation(1,1))^2+(intersection(1,2)-right_waterfall_deviation(1,2))^2);


 
 deviation_method_k1_and_k2 = [k1 k2];

% Luis Alberto Cañizares
