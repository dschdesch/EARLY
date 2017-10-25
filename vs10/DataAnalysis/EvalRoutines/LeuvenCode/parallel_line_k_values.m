% This method, finds the tangents of both sides of the "mountain" looking
% function and by drawing parallel lines to these points, it is possible to
% know where the waterfall shaped corner starts to form. The distance from
% these points to the tangents is known as the coefficients k1 and k2.
% Please note that the tangents are linear fits due to the surface
% roughness so an approximation of the coefficients k1 and k2 is only
% possible.

function parallel_lines_k1_and_k2 = parallel_line_k_values(ave_profile,length_profile, sliders)

global left_high
global right_high; 
global left_low;
global right_low;
global left_waterfall_parallel;
global right_waterfall_parallel;
global intersection;
global left_arm;
global right_arm;
global left_points;
global right_points;
global left_parallel;
global right_parallel;
global left_line;
global right_line;

global top_left_slider
global top_right_slider
global bottom_left_slider
global bottom_right_slider
              
top_left_slider = sliders(1,1);
top_right_slider = sliders(1,2);
bottom_left_slider = sliders(2,1);
bottom_right_slider = sliders(2,2);


% Find Max point
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
right_arm = ave_profile(i+right_high_clearance:length_profile-j2-left_low_clearance,1:2);


% Fitting lines through these points 
left_fit = polyfit(left_arm(:,1),left_arm(:,2),1);
left_points = ave_profile((j+left_low_clearance):(i+right_high_clearance),1);
left_line = polyval(left_fit,left_points);

right_fit = polyfit(right_arm(:,1),right_arm(:,2),1);
right_points = ave_profile((i-left_high_clearance):(length_profile - j2 - right_low_clearance),1);
right_line = polyval(right_fit,right_points);

% Parallel lines

global paralleltolerance
left_clearance = [0 paralleltolerance];                 % Parallel line distance
left_para_fit = left_fit(1,1:2)-left_clearance;         % Parallel line, ie. same slope diff y intersection point 
left_parallel = polyval(left_para_fit,left_points);     % Points for parallel line


right_clearance = left_clearance;                       % Same as above but 
right_para_fit = right_fit(1,1:2)-right_clearance;      % on the right side
right_parallel = polyval(right_para_fit,right_points);  %

% Find Waterfall Points
% These points can only be approximations due to surface roughness.
% LEFT
length_left = length(left_parallel);
count2 = i + right_high_clearance;

for count = length_left:-1:1
    if left_parallel(count,1) < ave_profile(count2, 2) 
        break
    end
 count2 = count2 - 1;   
end
left_waterfall_parallel = [ave_profile(count2,1) ave_profile(count2,2)];

% Right
length_right = length(right_parallel);
count3 = i - left_high_clearance;
for count = 1:length_right
    if right_parallel(count,1) < ave_profile(count3, 2) 
        break
    end
 count3 = count3 + 1;   
end
right_waterfall_parallel = [ave_profile(count3,1) ave_profile(count3,2)];


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
%


%
% Finding k1 and k2
%

k1 = sqrt((intersection(1)-left_waterfall_parallel(1))^2 +...
    (intersection(2)-left_waterfall_parallel(2))^2 - (left_clearance(2))^2);

k2 = sqrt((intersection(1)- right_waterfall_parallel(1))^2 +...
    (intersection(2)- right_waterfall_parallel(2))^2 -(right_clearance(2))^2);


parallel_lines_k1_and_k2 = [k1 k2];


% Luis Alberto Canizares