% This function is just for plotting a particular circle in the desired axes.


function plotleftcircle(handles)

global ave_profile;
global length_profile;
global intersection;


 axes(handles.axes3)
   [xfitl,yfitl,Rfitl,xfitr,yfitr,Rfitr,xfittotal,yfittotal,Rfittotal] = circlefit(ave_profile,length_profile);
           % LEFT Circle FIT    
           leftcircle = rectangle('position',[xfitl-Rfitl,yfitl-Rfitl,Rfitl*2,Rfitl*2],...
                   'curvature',[1,1],'linestyle','-','edgecolor','r');

                        plot(xfitl,yfitl,'g.')
                        xlim([xfitl-Rfitl-2,xfitl+Rfitl+2])
                        ylim([yfitl-Rfitl-2,yfitl+Rfitl+2])
           
          
               zoomcenter(intersection(1,1),intersection(1,2), 8)         
end

% Luis Alberto Canizares