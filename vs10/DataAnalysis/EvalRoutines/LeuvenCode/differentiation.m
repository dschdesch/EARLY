% This function uses a method of differentiation previously developed in TU
% Berlin. Once the differentiation has been done, the function finds and
% recognises different slopes caused by the differentiation of the curving 
% of the edge. 



function differentiation(Datasety,Datasetx,handles)
% Set up the Work Parameters
pA                  = 150;      %Amount of desired points
pA_spline           = 300;      %Amount of desired points for the spline interpolation
small_area_size    = round(length(Datasety)/pA);
x_step_size         = 0.001;

global ave_profile
global n

% Work in progress
lastx = 0;
sx(pA) = 0;
sy(pA) = 0;

for k=1:pA
    newx        = lastx + small_area_size*x_step_size-x_step_size;
    x           = lastx : 0.001 : newx;
    sx(k)       = sum(x) / length(x);
    beginData   = (k-1)*small_area_size + 1;
    endData     = (k)*small_area_size;
    sy(k)       = sum(Datasety(beginData : endData)) / length(Datasety(beginData : endData));
    lastx       = newx;
end
%Datasetx =(0:x_step_size:length(Datasety)*x_step_size - x_step_size)';

% Derivatives
for k=1:pA/2-1
derivative_y(k) =(sy((k+1)*2)-sy(k*2))/(sx((k+1)*2)-sx(k*2));
derivative_x(k) =(sx((k+1)*2)+sx(k*2))/2;        
end

% Get the points of interest (first the first, then last point, then all in between)
lengthderivativex = length(derivative_x);
countpoints = 0;
countstars = 0;


for count = 2:lengthderivativex
countpoints = countpoints + 1;
heightdifference = derivative_y(count-1) - derivative_y(count);

    if heightdifference > 0
     countstars = countstars + 1;
    end

    if heightdifference < 0
        countstars = 0;
    end

    if countstars >= 4
        break
    end


 end
countpoints = countpoints - 3;
topslopepoint = [derivative_x(countpoints),derivative_y(countpoints)];

% Counting backwards to get the last point
lengthderivativexb = length(derivative_x);   
countpointsb = 0;
countstarsb = 0;


for count = lengthderivativexb-1:-1:1
countpointsb = countpointsb + 1;
heightdifferenceb = derivative_y(count-1) - derivative_y(count);

    if heightdifferenceb > 0
     countstarsb = countstarsb + 1;
    end

    if heightdifferenceb < 0
        countstarsb = 0;
    end

    if countstarsb >= 4
        break
    end


end

lateraltolerance = 0.0005;
countpointsb = lengthderivativexb - countpointsb + 3;
checklateral = derivative_y(countpointsb-1) - derivative_y(countpointsb);

while checklateral < lateraltolerance       
% Check for lateral movements  
if checklateral < lateraltolerance
    countpointsb = countpointsb -1;
end
checklateral = derivative_y(countpointsb-1) - derivative_y(countpointsb);

end


bottomslopepoint = [derivative_x(countpointsb),derivative_y(countpointsb)];


% The rest of the points
skiingpoints(:,1) = derivative_x(countpoints:countpointsb);
skiingpoints(:,2) = derivative_y(countpoints:countpointsb);

totalslopepoints = skiingpoints;



% Finding the optimal amount of slopes.

maxnumberslopes = 6;                % This actually means 5 max slopes.
                                    % Program will automatically delete one slope.

lenghttotalslopepoints = length(totalslopepoints);
    if lenghttotalslopepoints <= maxnumberslopes
        n = lenghttotalslopepoints - 2;
    end

    if lenghttotalslopepoints > maxnumberslopes
        n = maxnumberslopes;
    end

nstick =  n;


% Fitting max amount of slopes
ab = BrokenStickRegression(totalslopepoints(:,1), totalslopepoints(:,2), nstick);
slopes(n,2) = 0;
for count = 1:n
   slopes(count,1:2) = polyfit(ab(count:count+1,1),ab(count:count+1,2),1);
end

% Finding tolerance for specific function
toleranceslopes(:,1) = abs(slopes(1:end-1,1) - slopes(2:end,1));
toleranceslopes = sum(toleranceslopes)/length(toleranceslopes);
for count = 2:n
    slopesdifference(count-1,1) = abs(slopes(count,1) - slopes(count-1,1));      
end
lengthslopesdifference = length(slopesdifference);
for count = 1:lengthslopesdifference
    if nstick <= 0
    nstick = 1;
    break
    end
    if slopesdifference(count,1) < toleranceslopes
        nstick = nstick - 1;
    end
end

% Dismissing the number of slopes that do not meet the tolerance
% requirements
nstick = nstick - 1;
if nstick <=1
    nstick = 1;
end


% Finding values of different slopes, this time with the correct amount of slopes.
ab = BrokenStickRegression(totalslopepoints(:,1), totalslopepoints(:,2), nstick);
for count = 1:nstick
   slopes(count,1:2) = polyfit(ab(count:count+1,1),ab(count:count+1,2),1);
end


% Plot
axes(handles.axes11)
cla(handles.axes11,'reset')
grid on;  hold on;

plot(derivative_x,derivative_y,'*-')

plot(skiingpoints(:,1),skiingpoints(:,2),'g*')
plot(bottomslopepoint(1,1),bottomslopepoint(1,2),'r*')
plot(topslopepoint(1,1),topslopepoint(1,2),'r*')
plot(ab(:,1),ab(:,2),'r-o')



axes(handles.axes12)
cla(handles.axes12,'reset')
grid on; axis equal;  hold on;
title('Average Profile');
xlabel('X');
ylabel('Z');
set(gca,'fontsize',11);
set(gcf,'color','white');
plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);


information ={' ' ' ';
              ' ' ' ';
              ' ' ' ';
              ' ' ' ';
              ' ' ' '; 
                     };
for count = 1:nstick
information(count, 1:2) = {sprintf('Slope %d',count) num2str(slopes(count,1))
    };
               
end

   
set(handles.uitable3,'Data', information)
set(handles.uitable3, 'ColumnWidth', {300 300});

% Luis Alberto Canizares (Except the derivative method)


