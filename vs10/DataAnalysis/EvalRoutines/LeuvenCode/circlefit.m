% This function finds the optimal amount of points for for then fitting a
% circle through them. 

% This gets the points and then it fits a circle trhough them 
function [xfitl,yfitl,Rfitl,xfitr,yfitr,Rfitr,xfittotal,yfittotal,Rfittotal] = circlefit(ave_profile,length_profile)

global intersection;
global left_waterfall_parallel;
left_pc = 0;
i = 0;
for count = 1:length_profile

    if ave_profile(count) >= left_waterfall_parallel
        i = i+1;   
        left_pc(i,1:2) = ave_profile(count,1:2);

    end
    
    if ave_profile(count,1) >= intersection(1,1)
       break
    end
end

global right_waterfall_parallel
right_pc = 0;
i = 0;
for count = 1:length_profile
    
    if ave_profile(count,1) >= intersection(1,1)
      i = i+1;   
        right_pc(i,1:2) = ave_profile(count,1:2);
    end
    
    if ave_profile(count) >= right_waterfall_parallel
        break
    end
    
end

[xfitl,yfitl,Rfitl] = circfit(left_pc(:,1),left_pc(:,2));
[xfitr,yfitr,Rfitr] = circfit(right_pc(:,1),right_pc(:,2));


% TOTAL Circle Fit

total_pc = 0;
i = 0; 

for count = 1:length_profile
    
    if ave_profile(count) >= left_waterfall_parallel
        i = i+1;   
        total_pc(i,1:2) = ave_profile(count,1:2);

    end
    
    if ave_profile(count,1) >=  right_waterfall_parallel
       break
    end
   
end

[xfittotal,yfittotal,Rfittotal] = circfit(total_pc(:,1),total_pc(:,2));

% Luis Alberto Canizares