% This function updates some of the data in the interface. This was done
% for speed and convinience. 


function updateinfo(sliders,handles) 
  
% Get Average Profile 
global ave_profile
global data_straight
global x;
global y;
global z;
global yRot;
global zRot;

ave_profile = average_profile(data_straight);

% [xx,yy] = testsingleprofile;
% ave_profile = [xx,yy];


% Choose a Profile. You can also chose a single profile instead of using
% an average one
% chosen = 200;                                     % Profile nº
% data_straight_2 = sortrows(data_straight,1);
% chosen_profile = data_straight_2(chosen*483:(chosen*482)+482,2:3);
% 
% ave_profile = chosen_profile;                   % For convinience remove
%                                                 % comments and add comments
%                                                 % above for a desired
                                                % profile

global length_profile
length_profile = length(ave_profile); 

          
% Get K1 and K2 from the parallel line method
parallel_lines_k1_and_k2 = parallel_line_k_values(ave_profile,length_profile,sliders);

% Get K1 and K2 from the deviation of line method
deviation_method_k1_and_k2 = deviation_method_k_values(ave_profile,length_profile,sliders);

set(handles.k1parallel,'String',num2str( parallel_lines_k1_and_k2(1,1)));
set(handles.k2parallel,'String',num2str( parallel_lines_k1_and_k2(1,2)));

set(handles.k1devi,'String',num2str( deviation_method_k1_and_k2(1,1)));
set(handles.k2devi,'String',num2str( deviation_method_k1_and_k2(1,2)));

set(handles.popsmall,'Value',2);
set(handles.popbig,'Value',4);
                 
                global left_high
                global right_high 
                global left_low
                global right_low
                global left_waterfall_parallel
                global right_waterfall_parallel
                global left_waterfall_deviation
                global right_waterfall_deviation
                global intersection
                global left_arm
                global right_arm
                global left_points
                global right_points
                global left_parallel
                global right_parallel
                global left_line
                global right_line
                 
              
                
                
               axes(handles.axes1);
               
               cla(handles.axes1,'reset')
                hold off;
                
                hold on;
                grid on; box on;axis equal;
                title('Parallel lines Method');
                set(gca,'fontsize',11);
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gcf,'color','white');

                plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);
                plot(left_high(1,1),left_high(1,2), 'go')
                plot(right_high(1,1),right_high(1,2), 'go')
                plot(left_low(1,1),left_low(1,2), 'go')
                plot(right_low(1,1),right_low(1,2), 'go')
                plot(left_arm(:,1),left_arm(:,2), 'g-')
                plot(right_arm(:,1),right_arm(:,2), 'g-')
                plot(left_points,left_line,'r-')
                plot(right_points,right_line,'r-')
                plot(left_points,left_parallel,'m-')
                plot(right_points,right_parallel,'m-')
                plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
                plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
                plot(intersection(1,1),intersection(1,2),'g+')
               
               
                axes(handles.axes3);
               cla(handles.axes3,'reset')
                hold on;
                axis equal
                grid on; box on;
                set(gca,'fontsize',11);
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gcf,'color','white');

                
     [xfitl,yfitl,Rfitl,xfitr,yfitr,Rfitr,xfittotal,yfittotal,Rfittotal] = circlefit(ave_profile,length_profile);

                
                
                plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);
                plot(left_high(1,1),left_high(1,2), 'go')
                plot(right_high(1,1),right_high(1,2), 'go')
                plot(left_low(1,1),left_low(1,2), 'go')
                plot(right_low(1,1),right_low(1,2), 'go')
                plot(left_points,left_line,'r-')
                plot(right_points,right_line,'r-')
                plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*') 
                plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
                plot(intersection(1,1),intersection(1,2),'g+')
                
               
                
                
               xyfitl = [xfitl*1000,yfitl*1000];     % 1000 is to scale mm to microns
               xyfitr = [xfitr*1000,yfitr*1000];
               xyfittotal = [xfittotal*1000,yfittotal*1000];
               information =  {' Intersection Point' num2str(intersection);
                               ' Left Waterfall Parallel Method' num2str(left_waterfall_parallel);
                               ' Right Waterfall Parallel Method' num2str(right_waterfall_parallel);
                               ' Left Waterfall Deviation Method' num2str(left_waterfall_deviation);
                               ' Right Waterfall Deviation Method' num2str(right_waterfall_deviation);
                               ' ' ' ';
                               ' RADIUS OF CIRCLE'  '(microns)';
                               ' Left Circle' num2str(Rfitl*1000);
                               ' Right Circle' num2str(Rfitr*1000);
                               ' Total Circle' num2str(Rfittotal*1000);
                               ' ' ' ';
                               ' CENTER OF CIRCLE' 'Coordinates (x,z) (microns)';
                               ' Left Circle' num2str(xyfitl);
                               ' Right Circle' num2str(xyfitr);
                               ' Total Circle' num2str(xyfittotal);
                               ' ' ' ';
                               ' ' ' ';
                               ' K1 & K2' '(microns)';
                               ' PARALLEL METHOD' ' ';
                               ' K1' num2str( parallel_lines_k1_and_k2(1,1)*1000);
                               ' K2' num2str( parallel_lines_k1_and_k2(1,2)*1000);
                               ' ' ' ';
                               ' DEVIATION METHOD ' ' ';
                               ' K1' num2str( deviation_method_k1_and_k2(1,1)*1000);
                               ' K2' num2str( deviation_method_k1_and_k2(1,2)*1000);
                               ' ' ' ';
                               ' ' ' ';
                               ' ' ' ';
                               ' ' ' ';
                               ' ' ' ';
                               ' ' ' '};
               
               
            
               set(handles.informationtable,'Data', information)
               set(handles.informationtable,'ColumnName', {'Point','Coordinates (x,z)'})
               set(handles.informationtable, 'ColumnWidth', {220 2200});
              
               
             axes(handles.axes4);
               cla(handles.axes4,'reset')
               grid on; box on;
                hold on;   
               axis equal
                plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);
                plot(left_high(1,1),left_high(1,2), 'go')
                plot(right_high(1,1),right_high(1,2), 'go')
                plot(left_low(1,1),left_low(1,2), 'go')
                plot(right_low(1,1),right_low(1,2), 'go')
                plot(left_points,left_line,'r-')
                plot(right_points,right_line,'r-')
         
end
    

% Luis Alberto Canizares