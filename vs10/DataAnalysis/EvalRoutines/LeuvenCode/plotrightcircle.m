% This function is just for plotting a particular circle in the desired axes. 

function plotrightcircle(handles)
global ave_profile;
global length_profile;
global intersection

 axes(handles.axes3)
   [xfitl,yfitl,Rfitl,xfitr,yfitr,Rfitr,xfittotal,yfittotal,Rfittotal] = circlefit(ave_profile,length_profile);
          
            % RIGHT Circle FIT            
            rightcircle =  rectangle('position',[xfitr-Rfitr,yfitr-Rfitr,Rfitr*2,Rfitr*2],...
                   'curvature',[1,1],'linestyle','-','edgecolor','r');

                        plot(xfitr,yfitr,'g.')
                        xlim([xfitr-Rfitr-2,xfitr+Rfitr+2])
                        ylim([yfitr-Rfitr-2,yfitr+Rfitr+2])
          
       zoomcenter(intersection(1,1),intersection(1,2), 8)                 

end

% Luis Alberto Canizares 