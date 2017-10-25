% This function is just for plotting a particular circle in the desired axes. 

function plottotalcircle(handles)

global ave_profile;
global length_profile;
global intersection

 axes(handles.axes3)
   [xfitl,yfitl,Rfitl,xfitr,yfitr,Rfitr,xfittotal,yfittotal,Rfittotal] = circlefit(ave_profile,length_profile);
                  
            % TOTAL Circle FIT        
             totalcircle = rectangle('position',[xfittotal-Rfittotal,yfittotal - Rfittotal,Rfittotal*2,Rfittotal*2],...
                   'curvature',[1,1],'linestyle','-','edgecolor','g');

                        plot(xfittotal,yfittotal,'g.')
                        xlim([xfittotal-Rfittotal-2,xfittotal+Rfittotal+2])
                        ylim([yfittotal-Rfittotal-2,yfittotal+Rfittotal+2])
                        
   zoomcenter(intersection(1,1),intersection(1,2), 8)
                        
end

% Luis Alberto Canizares