% Luis Alberto Canizares
% This m-file together with Interface.fig contains the interface body and
% design. This means that this is the core of the program. Any other m-file
% will have to be called using these interface's figure callback's.



function varargout = Interface(varargin)
% INTERFACE MATLAB code for Interface.fig
%      INTERFACE, by itself, creates a new INTERFACE or raises the existing
%      singleton*.
%
%      H = INTERFACE returns the handle to a new INTERFACE or the handle to
%      the existing singleton*.
%
%      INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACE.M with the given input arguments.
%
%      INTERFACE('Property','Value',...) creates a new INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Interface

% Last Modified by GUIDE v2.5 08-Aug-2012 10:00:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Interface_OpeningFcn, ...
                   'gui_OutputFcn',  @Interface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT
end

% --- Executes just before Interface is made visible.
function Interface_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for Interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
%This function runs just after the program starts running. ie. all initial
%configurations must go here.
function varargout = Interface_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


    set(handles.controlpanel,'Visible', 'on')
    set(handles.startpanel, 'Visible', 'on')
    set(handles.circlefitpanel,'Visible', 'off')
    set(handles.calibrationpanel,'Visible','off')
    set(handles.differentiationpanel,'Visible','off')
    set(handles.calibrationbutton,'Visible', 'off')      % Not recomended to use unless is strictly necessary

    set(handles.text49,'Visible','off')
    set(handles.text50,'Visible','off')
    set(handles.text51,'Visible','off')
    set(handles.text52,'Visible','off')
    set(handles.text53,'Visible','off')
    set(handles.text54,'Visible','off')
    set(handles.text55,'Visible','off')
    set(handles.text56,'Visible','off')
    set(handles.text57,'Visible','off')
    set(handles.text58,'Visible','off')
    set(handles.text59,'Visible','off')
    set(handles.text60,'Visible','off')
    set(handles.text61,'Visible','off')
    set(handles.text62,'Visible','off') 
    
    set(handles.text63,'Visible','on')
    set(handles.text64,'Visible','off')
    set(handles.text65,'Visible','off')
    set(handles.text66,'Visible','off') 
    



information1 = {' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' '; 
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' ';
    ' ' ' '};
               set(handles.informationtable,'Data', information1)
               set(handles.informationtable,'ColumnName', {'Point','Coordinates (x,z)'})
               set(handles.informationtable, 'ColumnWidth', {200 200});



information2 ={' ' ' ';
              ' ' ' ';
              ' ' ' ';
              ' ' ' ';
              ' ' ' '; 
                     };
set(handles.uitable3,'Data', information2)
set(handles.uitable3, 'ColumnWidth', {250 250});

end


% --- Executes on button press in load_data.

% This function loads the data and then updates the information.
function load_data_Callback(hObject, eventdata, handles)

[FileName,PathName] = uigetfile('*.txt','SelectFiles','MultiSelect','Off');                            
% Checks done in this order!!
if ~iscell(FileName) && (length(FileName) == 1) && (FileName == 0)       % User canceled
    return;
else
    combinedStr = strcat(PathName,FileName); 
    fid = fopen(combinedStr);                                            % Name of the txt file
    handles.data = textscan(fid, '%f %f %f','headerLines', 9);           % Last number indicates how 
    fclose('all');                                                       % many lines will matlab skip
                                                                         % from the txt file. ie, headerlines.
end

set(handles.text136,'String', FileName)    

global x;
global y;
global z;
global yRot;
global zRot;

x = cell2mat(handles.data(:,1));
y = cell2mat(handles.data(:,2));
z = cell2mat(handles.data(:,3));
 

 [x,yRot,zRot] = straight_points(x, y, z);         % This function organizes 
                                                   % the data so it is aligned to the x,y 
                                                   % and z axes, in order to minimise errors
 
 global data_straight
 data_straight = [x yRot zRot];
 
 global top_left_slider
 global top_right_slider
 global bottom_left_slider
 global bottom_right_slider

   top_left_slider = 40;
   top_right_slider = 50;
   bottom_left_slider = 60;
   bottom_right_slider = 60;
  
  
 global sliders
   sliders = [top_left_slider      top_right_slider;
              bottom_left_slider   bottom_right_slider];
global paralleltolerance
global devitolerance

paralleltolerance = 0.001;
devitolerance = 0.001;


   updateinfo(sliders,handles);                 
   
                 axes(handles.axes2);
                 
                 tri = delaunay(x,yRot);
                  trisurf(tri, x,yRot,zRot);
                  title('Profile Aligned')
                  xlabel('x'); ylabel('y'); zlabel('z');
                  axis equal;
                  set(gca,'fontsize',11);
                  shading interp
                  caxis([0 0.6]), colormap('default'), colorbar
    
                  
    set(handles.topleftreference,'Min',-50,'Max',0,...
        'SliderStep',[0.01 0.1],'Value',-20);
    set(handles.toprightreference,'Min',-50,'Max',0,...
        'SliderStep',[0.01 0.1],'Value',-20);
    set(handles.bottomleftreference,'Min',0,'Max',100,...
        'SliderStep',[0.01 0.1],'Value',20);
    set(handles.bottomrightreference,'Min',0,'Max',100,...
        'SliderStep',[0.01 0.1],'Value',20);
    
    set(handles.texttl,'String', num2str(sliders(1,1)))
    set(handles.textbl,'String', num2str(sliders(2,1)))
    set(handles.texttr,'String', num2str(sliders(1,2)))
    set(handles.textbr,'String', num2str(sliders(2,2)))
    
    
    set(handles.deviationtoleranceslider, 'Min', 0,'Max', 0.02,...
        'SliderStep',[0.01 0.1],'Value', 0.001);
    set(handles.paralleltoleranceslider, 'Min', 0,'Max', 0.02,...
        'SliderStep',[0.01 0.1],'Value', 0.001);
    
    set(handles.edittolerancedevi, 'String',...
        num2str(get(handles.deviationtoleranceslider,'Value')))
    set(handles.edittolerancepara, 'String',...
        num2str(get(handles.deviationtoleranceslider,'Value')))
    
    
   differentiation(zRot,x,handles);
    
% Update handles structure
guidata(hObject, handles);
end
% --- Executes on selection change in popbig.
function popbig_Callback(hObject, eventdata, handles)

     axes(handles.axes1);
        str = get(hObject, 'String');
        val = get(hObject,'Value');
        
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
                global ave_profile
                global x
                global y
                global z
                global yRot
                global zRot
         % Set current data to the selected data set.
         switch str{val};
             case '3D Original'                                 % User selects 3D Original
                
                 hold off;  
                 tri = delaunay(x,y);
                  trisurf(tri, x,y,z);
                  title('Profile Original');
                  xlabel('x'); ylabel('y'); zlabel('z');
                  set(gca,'fontsize',11);
                  axis equal;
                  shading interp;
                  caxis([-0.3 0.3]), colormap('default'), colorbar;
                hold off;                
                
             case '3D Aligned'                                   % User selects 3D aligned
               
                 hold off;
                 tri = delaunay(x,yRot);
                  trisurf(tri, x,yRot,zRot);
                  title('Profile Aligned')
                  xlabel('x'); ylabel('y'); zlabel('z');
                  axis equal;
                  set(gca,'fontsize',11);
                  shading interp
                  caxis([0 0.6]), colormap('default'), colorbar
               hold off;
            
             case 'Average Profile'                              % User selects Average Profile
               cla(handles.axes1,'reset')
                grid on; axis equal;  hold on;
                title('Average Profile');
                xlabel('X');
                ylabel('Z');
                set(gca,'fontsize',11);
                set(gcf,'color','white');
                plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);
               hold off;

             case 'Parallel Lines Method'                        % User selects the method of Parallel lines
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
  
             case 'Deviation Method'                             % User selects the deviation method

                 cla(handles.axes1,'reset')
                 hold off;
                 
                   hold on; grid on; box on;
                    title('Deviation Method');
                    xlabel('X (mm)');
                    ylabel('Z (mm)');
                    set(gca,'fontsize',11);
                    set(gcf,'color','white');
                    axis equal;
                    plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);
                    plot(left_high(1,1),left_high(1,2), 'go')
                    plot(right_high(1,1),right_high(1,2), 'go')
                    plot(left_low(1,1),left_low(1,2), 'go')
                    plot(right_low(1,1),right_low(1,2), 'go')
                    plot(left_arm(:,1),left_arm(:,2), 'g-')
                    plot(right_arm(:,1),right_arm(:,2), 'g-')
                    plot(left_points(:,1),left_line,'r-')
                    plot(right_points(:,1),right_line,'r-')
                    plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*')
                    plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*')
                    plot(intersection(1,1),intersection(1,2),'g+')


         end


end

% --- Executes during object creation, after setting all properties.
function popbig_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popbig.
function popsmall_Callback(hObject, eventdata, handles)

axes(handles.axes2);
str = get(hObject, 'String');
       val = get(hObject,'Value');
        
                
                
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
                global ave_profile
                global x
                global y
                global z
                global yRot
                global zRot
         % Set current data to the selected data set.
         switch str{val};
             case '3D Original'                                 % User selects 3D Original
                
                 hold off;
                 tri = delaunay(x,y);
                  trisurf(tri, x,y,z);
                  title('Profile Original')
                  xlabel('x'); ylabel('y'); zlabel('z');
                  set(gca,'fontsize',11);
                  axis equal;
                  shading interp
                  caxis([-0.3 0.3]), colormap('default'), colorbar
                hold off;
                
             case '3D Aligned'                                   % User selects 3D aligned
               
                 hold off;
                 tri = delaunay(x,yRot);
                  trisurf(tri, x,yRot,zRot);
                  title('Profile Aligned')
                  xlabel('x'); ylabel('y'); zlabel('z');
                  set(gca,'fontsize',11);
                  axis equal;
                  shading interp
                  caxis([0 0.6]), colormap('default'), colorbar
               hold off;
            
             case 'Average Profile'                              % User selects Average Profile
               cla(handles.axes2,'reset')
                grid on; axis equal;  hold on;
                title('Average Profile');
                xlabel('X');
                ylabel('Z');
                set(gca,'fontsize',11);
                set(gcf,'color','white');
                plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);
               hold off;

             case 'Parallel Lines Method'                        % User selects the method of Parallel lines
              cla(handles.axes2,'reset')
                hold off;
                
                hold on;
                grid on; box on;axis equal;
                title('Parallel lines Method');
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gca,'fontsize',11);
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
  
             case 'Deviation Method'                             % User selects the deviation method

                 cla(handles.axes2,'reset')
                 hold off;
                 
                   hold on; grid on; box on;
                    title('Deviation Method');
                    xlabel('X (mm)');
                    ylabel('Z (mm)');
                    set(gca,'fontsize',11);
                    set(gcf,'color','white');
                    axis equal;
                    plot(ave_profile(:,1), ave_profile(:,2), '-','linewidth',2);
                    plot(left_high(1,1),left_high(1,2), 'go')
                    plot(right_high(1,1),right_high(1,2), 'go')
                    plot(left_low(1,1),left_low(1,2), 'go')
                    plot(right_low(1,1),right_low(1,2), 'go')
                    plot(left_arm(:,1),left_arm(:,2), 'g-')
                    plot(right_arm(:,1),right_arm(:,2), 'g-')
                    plot(left_points(:,1),left_line,'r-')
                    plot(right_points(:,1),right_line,'r-')
                    plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*')
                    plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*')
                    plot(intersection(1,1),intersection(1,2),'g+')


         end

end
         
% --- Executes during object creation, after setting all properties.
function popsmall_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
    end
end

% --- Executes during object creation, after setting all properties.
function k1parallel_CreateFcn(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
set(handles.startpanel,'Visible','on')
set(handles.circlefitpanel,'Visible','off')
set(handles.calibrationpanel,'Visible','off')
set(handles.differentiationpanel,'Visible','off')

end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
set(handles.startpanel,'Visible','off')
set(handles.circlefitpanel,'Visible','on')
set(handles.calibrationpanel, 'Visible','off')
set(handles.differentiationpanel,'Visible','off')
    set(handles.text63,'Visible','off')
    set(handles.text64,'Visible','on')
    set(handles.text65,'Visible','off')
    set(handles.text66,'Visible','off')
     
   
end


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function informationtable_CreateFcn(hObject, eventdata, handles)
end


% --- Executes when selected cell(s) is changed in informationtable.
function informationtable_CellSelectionCallback(hObject, eventdata, handles)
end


% --------------------------------------------------------------------
function uipushtool8_ClickedCallback(hObject, eventdata, handles)
load_data_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in calibrationbutton.
function calibrationbutton_Callback(hObject, eventdata, handles)
set(handles.startpanel,'Visible','off')
set(handles.circlefitpanel,'Visible','off')
set(handles.calibrationpanel, 'Visible','on')
set(handles.differentiationpanel,'Visible','off')

    set(handles.text63,'Visible','off')
    set(handles.text64,'Visible','off')
    set(handles.text65,'Visible','off')
    set(handles.text66,'Visible','on')


end


% --- Executes on slider movement.
function topleftreference_Callback(hObject, eventdata, handles)

global top_left_slider
slider_value = get(hObject,'Value');
slider_value = round(slider_value);
top_left_slider = -slider_value;

 global sliders
 sliders(1,1) = top_left_slider;
    updateinfo(sliders,handles);
    
 set(handles.texttl,'String', num2str(-slider_value))
 

            global left_waterfall_parallel
            global right_waterfall_parallel
            global intersection
            global left_points
            global right_points
            global left_parallel
            global right_parallel
            global left_waterfall_deviation
            global right_waterfall_deviation
 
 if get(handles.showparallel,'Value')  == get(handles.showparallel, 'Max')
            plot(left_points,left_parallel,'m-')
            plot(right_points,right_parallel,'m-')
            plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
            plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
            plot(intersection(1,1),intersection(1,2),'g+')

 end
 if get(handles.showdeviation,'Value')  == get(handles.showdeviation, 'Max')
            plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
            plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
            plot(intersection(1,1),intersection(1,2),'g+');
 end 
 
end

% --- Executes during object creation, after setting all properties.
function topleftreference_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function bottomleftreference_Callback(hObject, eventdata, handles)
global bottom_left_slider
slider_value = get(hObject,'Value');
slider_value = round(slider_value);
bottom_left_slider = slider_value;

 global sliders
 sliders(2,1) = bottom_left_slider;
    updateinfo(sliders,handles);
    
    set(handles.textbl,'String', num2str(slider_value))
            global left_waterfall_parallel
            global right_waterfall_parallel
            global intersection
            global left_points
            global right_points
            global left_parallel
            global right_parallel
            global left_waterfall_deviation
            global right_waterfall_deviation
    if get(handles.showparallel,'Value')  == get(handles.showparallel, 'Max')
            plot(left_points,left_parallel,'m-')
            plot(right_points,right_parallel,'m-')
            plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
            plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
            plot(intersection(1,1),intersection(1,2),'g+')

    end
 if get(handles.showdeviation,'Value')  == get(handles.showdeviation, 'Max')
            plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
            plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
            plot(intersection(1,1),intersection(1,2),'g+');
 end 
 
end

% --- Executes during object creation, after setting all properties.
function bottomleftreference_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function bottomrightreference_Callback(hObject, eventdata, handles)
global bottom_right_slider
slider_value = get(hObject,'Value');
slider_value = round(slider_value);
bottom_right_slider = slider_value;

 global sliders
 sliders(2,2) = bottom_right_slider;
    updateinfo(sliders,handles);
    
    set(handles.textbr,'String', num2str(slider_value))   
            global left_waterfall_parallel
            global right_waterfall_parallel
            global intersection
            global left_points
            global right_points
            global left_parallel
            global right_parallel
            global left_waterfall_deviation
            global right_waterfall_deviation
    
    if get(handles.showparallel,'Value')  == get(handles.showparallel, 'Max')
            plot(left_points,left_parallel,'m-')
            plot(right_points,right_parallel,'m-')
            plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
            plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
            plot(intersection(1,1),intersection(1,2),'g+')

    end
 if get(handles.showdeviation,'Value')  == get(handles.showdeviation, 'Max')
            plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
            plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
            plot(intersection(1,1),intersection(1,2),'g+');
 end 
 
end

% --- Executes during object creation, after setting all properties.
function bottomrightreference_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function toprightreference_Callback(hObject, eventdata, handles)
global top_right_slider
slider_value = get(hObject,'Value');
slider_value = round(slider_value);
top_right_slider = -slider_value;

 global sliders
 sliders(1,2) = top_right_slider;
    updateinfo(sliders,handles);
    set(handles.texttr,'String', num2str(-slider_value))   
            global left_waterfall_parallel
            global right_waterfall_parallel
            global intersection
            global left_points
            global right_points
            global left_parallel
            global right_parallel
            global left_waterfall_deviation
            global right_waterfall_deviation
    if get(handles.showparallel,'Value')  == get(handles.showparallel, 'Max')
            plot(left_points,left_parallel,'m-')
            plot(right_points,right_parallel,'m-')
            plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
            plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
            plot(intersection(1,1),intersection(1,2),'g+')

    end
 if get(handles.showdeviation,'Value') == get(handles.showdeviation, 'Max')
            plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
            plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
            plot(intersection(1,1),intersection(1,2),'g+');
 end 
 
end

% --- Executes during object creation, after setting all properties.
function toprightreference_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on button press in showparallel.
function showparallel_Callback(hObject, eventdata, handles)
            
            global ave_profile
            global left_high
            global right_high 
            global left_low
            global right_low
            global left_waterfall_parallel
            global right_waterfall_parallel
            global intersection
            global left_points
            global right_points
            global left_line
            global right_line
            global left_parallel
            global right_parallel
            global left_waterfall_deviation
            global right_waterfall_deviation
axes(handles.axes4)
grid on; box on;
if get(hObject,'Value') == get(hObject, 'Max')
  
        plot(left_points,left_parallel,'m-')
        plot(right_points,right_parallel,'m-')
        plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
        plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
        plot(intersection(1,1),intersection(1,2),'g+')

end

if get(hObject,'Value') == get(hObject, 'Min')
        
       axes(handles.axes4);
      
               cla(handles.axes4,'reset')
                
               hold on;   
                grid on; box on;
               axis equal
                plot(ave_profile(:,1),ave_profile(:,2), '-','linewidth',2);
                plot(left_high(1,1),left_high(1,2), 'go')
                plot(right_high(1,1),right_high(1,2), 'go')
                plot(left_low(1,1),left_low(1,2), 'go')
                plot(right_low(1,1),right_low(1,2), 'go')
                plot(left_points,left_line,'r-')
                plot(right_points,right_line,'r-')
         if get(handles.showdeviation,'Value')  == get(handles.showdeviation, 'Max')
            plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
            plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
            plot(intersection(1,1),intersection(1,2),'g+');
         end 

end

end
% --- Executes on button press in showdeviation.
function showdeviation_Callback(hObject, eventdata, handles)
                global ave_profile
                global left_high
                global right_high 
                global left_low
                global right_low
                global left_waterfall_deviation
                global right_waterfall_deviation
                global intersection
                global left_points
                global right_points
                global left_line
                global right_line
                global left_waterfall_parallel
                global right_waterfall_parallel
                global left_parallel
                 global right_parallel

axes(handles.axes4)
grid on; box on;
if get(hObject,'Value') == get(hObject, 'Max')
   plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
   plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
   plot(intersection(1,1),intersection(1,2),'g+');
end

if get(hObject,'Value') == get(hObject, 'Min')
        
       axes(handles.axes4);
               cla(handles.axes4,'reset')
               grid on; box on;
               hold on;   
               axis equal
                plot(ave_profile(:,1),ave_profile(:,2), '-','linewidth',2);
                plot(left_high(1,1),left_high(1,2), 'go')
                plot(right_high(1,1),right_high(1,2), 'go')
                plot(left_low(1,1),left_low(1,2), 'go')
                plot(right_low(1,1),right_low(1,2), 'go')
                plot(left_points,left_line,'r-')
                plot(right_points,right_line,'r-')
         if get(handles.showparallel,'Value')  == get(handles.showparallel, 'Max')
            plot(left_points,left_parallel,'m-')
            plot(right_points,right_parallel,'m-')
            plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
            plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
            plot(intersection(1,1),intersection(1,2),'g+')

         end
 end

end


% --- Executes on button press in helpcalibration.
function helpcalibration_Callback(hObject, eventdata, handles)
end

% --- Executes on slider movement.
function paralleltoleranceslider_Callback(hObject, eventdata, handles)
global sliders
global paralleltolerance
slider_value = get(hObject,'Value');

paralleltolerance = slider_value;

    updateinfo(sliders,handles);
    set(handles.edittolerancepara,'String', num2str(slider_value))   
            
            global left_waterfall_parallel
            global right_waterfall_parallel
            global intersection
            global left_points
            global right_points
            global left_parallel
            global right_parallel
            global left_waterfall_deviation
            global right_waterfall_deviation
    if get(handles.showparallel,'Value')  == get(handles.showparallel, 'Max')
            plot(left_points,left_parallel,'m-')
            plot(right_points,right_parallel,'m-')
            plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
            plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
            plot(intersection(1,1),intersection(1,2),'g+')

    end
 if get(handles.showdeviation,'Value') == get(handles.showdeviation, 'Max')
            plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
            plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
            plot(intersection(1,1),intersection(1,2),'g+');
 end 
end

% --- Executes during object creation, after setting all properties.
function paralleltoleranceslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end

% --- Executes on slider movement.
function deviationtoleranceslider_Callback(hObject, eventdata, handles)

global devitolerance
slider_value = get(hObject,'Value');

devitolerance = slider_value;
global sliders
    updateinfo(sliders,handles);
    set(handles.edittolerancedevi,'String', num2str(slider_value))   
            
            global left_waterfall_parallel
            global right_waterfall_parallel
            global intersection
            global left_points
            global right_points
            global left_parallel
            global right_parallel
            global left_waterfall_deviation
            global right_waterfall_deviation
    if get(handles.showparallel,'Value')  == get(handles.showparallel, 'Max')
            plot(left_points,left_parallel,'m-')
            plot(right_points,right_parallel,'m-')
            plot(left_waterfall_parallel(1,1),left_waterfall_parallel(1,2),'r*')
            plot(right_waterfall_parallel(1,1),right_waterfall_parallel(1,2),'r*')
            plot(intersection(1,1),intersection(1,2),'g+')

    end
 if get(handles.showdeviation,'Value') == get(handles.showdeviation, 'Max')
            plot(left_waterfall_deviation(1,1),left_waterfall_deviation(1,2), 'r*');
            plot(right_waterfall_deviation(1,1),right_waterfall_deviation(1,2),'r*');
            plot(intersection(1,1),intersection(1,2),'g+');
 end 


end

% --- Executes during object creation, after setting all properties.
function deviationtoleranceslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in totalcirctoggle.
function totalcirctoggle_Callback(hObject, eventdata, handles)
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
global ave_profile

if get(hObject,'Value') == get(hObject,'Max')
    
    plottotalcircle(handles);
 
end
if get(hObject,'Value') == get(hObject,'Min')
    axes(handles.axes3);
               cla(handles.axes3,'reset')
                hold on;
                axis equal
                grid on; box on;
                set(gca,'fontsize',11);
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gcf,'color','white');

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

                
                if get(handles.rightcirctoggle, 'Value') == get(handles.rightcirctoggle,'Max')
                    plotrightcircle(handles)
                end
                if get(handles.leftcirctoggle, 'Value') == get(handles.leftcirctoggle,'Max')
                    plotleftcircle(handles)
                end
                
end

 if get(handles.rightcirctoggle, 'Value') == get(handles.rightcirctoggle,'Max') && ...
         get(handles.leftcirctoggle, 'Value') == get(handles.leftcirctoggle,'Max') && ...
         get(handles.totalcirctoggle, 'Value') == get(handles.totalcirctoggle,'Max')
                    
     zoomcenter(intersection(1,1),intersection(1,2), 3)         

 end
end
% --- Executes on button press in rightcirctoggle.
function rightcirctoggle_Callback(hObject, eventdata, handles)
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
global ave_profile

if get(hObject,'Value') == get(hObject,'Max')
    
    plotrightcircle(handles);
end
if get(hObject,'Value') == get(hObject,'Min')
    axes(handles.axes3);
               cla(handles.axes3,'reset')
                hold on;
                axis equal
                grid on; box on;
                set(gca,'fontsize',11);
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gcf,'color','white');

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
                
                if get(handles.totalcirctoggle, 'Value') == get(handles.totalcirctoggle,'Max')
                    plottotalcircle(handles)
                end
                if get(handles.leftcirctoggle, 'Value') == get(handles.leftcirctoggle,'Max')
                    plotleftcircle(handles)
                end
                
end



if get(handles.rightcirctoggle, 'Value') == get(handles.rightcirctoggle,'Max') && ...
         get(handles.leftcirctoggle, 'Value') == get(handles.leftcirctoggle,'Max') && ...
         get(handles.totalcirctoggle, 'Value') == get(handles.totalcirctoggle,'Max')
                    
     zoomcenter(intersection(1,1),intersection(1,2), 3)         

end
end
% --- Executes on button press in leftcirctoggle.
function leftcirctoggle_Callback(hObject, eventdata, handles)

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
global ave_profile

if get(hObject,'Value') == get(hObject,'Max')
    
    plotleftcircle(handles);
end
if get(hObject,'Value') == get(hObject,'Min')
    axes(handles.axes3);          
               cla(handles.axes3,'reset')
                hold on;
                axis equal
                grid on; box on;
                set(gca,'fontsize',11);
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gcf,'color','white');

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
                      
                if get(handles.rightcirctoggle, 'Value') == get(handles.rightcirctoggle,'Max')
                    plotrightcircle(handles)
                end
                if get(handles.totalcirctoggle, 'Value') == get(handles.totalcirctoggle,'Max')
                    plottotalcircle(handles)
                end
                
end

if get(handles.rightcirctoggle, 'Value') == get(handles.rightcirctoggle,'Min') && ...
         get(handles.leftcirctoggle, 'Value') == get(handles.leftcirctoggle,'Min') && ...
         get(handles.totalcirctoggle, 'Value') == get(handles.totalcirctoggle,'Min')
                  

end
end


% --- Executes on button press in differentiation.
function differentiation_Callback(hObject, eventdata, handles)

set(handles.controlpanel,'Visible', 'on')
set(handles.startpanel, 'Visible', 'off')
set(handles.circlefitpanel,'Visible', 'off')
set(handles.calibrationpanel,'Visible','off')
set(handles.differentiationpanel,'Visible','on')

    set(handles.text63,'Visible','off')
    set(handles.text64,'Visible','off')
    set(handles.text65,'Visible','on')
    set(handles.text66,'Visible','off')
    
   
end


% --- Executes on button press in togglehelp.
function togglehelp_Callback(hObject, eventdata, handles)
helpvalue = get(hObject, 'Value');


if helpvalue == 1;
    set(handles.text49,'Visible','on')
    set(handles.text50,'Visible','on')    
    set(handles.text51,'Visible','on')   
    set(handles.text52,'Visible','on')    
    set(handles.text53,'Visible','on')    
    set(handles.text54,'Visible','on')    
    set(handles.text55,'Visible','on')    
    set(handles.text56,'Visible','on')      
    set(handles.text57,'Visible','on')    
    set(handles.text58,'Visible','on')       
    set(handles.text59,'Visible','on')   
    set(handles.text60,'Visible','on')    
    set(handles.text61,'Visible','on')    
    set(handles.text62,'Visible','on')   
    
    
end

if helpvalue == 0;
    set(handles.text49,'Visible','off')    
    set(handles.text50,'Visible','off')    
    set(handles.text51,'Visible','off')    
    set(handles.text52,'Visible','off')    
    set(handles.text53,'Visible','off')    
    set(handles.text54,'Visible','off')    
    set(handles.text55,'Visible','off')    
    set(handles.text56,'Visible','off')
    set(handles.text57,'Visible','off')    
    set(handles.text58,'Visible','off')
    set(handles.text59,'Visible','off')    
    set(handles.text60,'Visible','off')    
    set(handles.text61,'Visible','off')    
    set(handles.text62,'Visible','off')
        
    
end
end


% --- Executes on button press in activatecalibration.
function activatecalibration_Callback(hObject, eventdata, handles)

calibrationvalue = get(hObject, 'Value');
if calibrationvalue == 0;
    set(handles.calibrationbutton,'Visible','off')
end
if calibrationvalue == 1;
    set(handles.calibrationbutton,'Visible','on')
end
end



function edittolerancedevi_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edittolerancedevi_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edittolerancepara_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edittolerancepara_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
