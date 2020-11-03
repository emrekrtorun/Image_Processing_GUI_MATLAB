function varargout = project(varargin)
% PROJECT MATLAB code for project.fig
%      PROJECT, by itself, creates a new PROJECT or raises the existing
%      singleton*.
%
%      H = PROJECT returns the handle to a new PROJECT or the handle to
%      the existing singleton*.
%
%      PROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT.M with the given input arguments.
%
%      PROJECT('Property','Value',...) creates a new PROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before project_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to project_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help project

% Last Modified by GUIDE v2.5 22-Jan-2019 13:10:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @project_OpeningFcn, ...
                   'gui_OutputFcn',  @project_OutputFcn, ...
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


% --- Executes just before project is made visible.
function project_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to project (see VARARGIN)

% Choose default command line output for project
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% UIWAIT makes project wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_image_Callback(hObject, eventdata, handles)

global im im2 baseFileName

[baseFileName,folder]=uigetfile('*') ;% read the input image
im=imread([folder,baseFileName]);

im2=im; %for backup process :)
axes(handles.axes1);
imshow(im);



% --------------------------------------------------------------------
function image_info_Callback(hObject, eventdata, handles)

global baseFileName

inf=imfinfo(baseFileName);
details= sprintf('Name = %s\n Width= %d\nHeight=%d \nColorType=%s\nBitsPerSample=%d\n',baseFileName, inf.Width, inf.Height, inf.ColorType, inf.NumberOfSamples);
set(handles.static,'string',details);

% --------------------------------------------------------------------
function save_image_Callback(hObject, eventdata, handles)

F = getframe(handles.axes2);
Image = frame2im(F);
% imwrite(Image,'C:/Users/Emre Can/Desktop/Matlab and Simulink/Project/output.jpg','jpg');
[F,P]= uiputfile();
imwrite(Image,fullfile(P,F));


% --- Executes on button press in plot_histogram.
function plot_histogram_Callback(hObject, eventdata, handles)

global im2;
image=im2;
[m,n,d]=size(image);
if (d==1)
    Channel1 = image(:,:,1);
    [M,N]=size(Channel1);
    t=1:256;
    n=0:255;
    count=0;
    for z=1:256
        for i=1:M
            for j=1:N
                if Channel1(i,j)==z-1
                    count=count+1;
                end
            end
        end
                t(z)=count;
                count=0;
    end
axes(handles.axes2);
plot(n,t,'Blue');
xlabel('Pixel Intensity');
ylabel('Number of Pixels');
xlim([0 256])

else
Red = image(:,:,1);
Green = image(:,:,2);
Blue = image(:,:,3);

[M,N]=size(Red);

t=1:256;
n=0:255;
count=0;
for z=1:256
    for i=1:M
        for j=1:N
            if Red(i,j)==z-1
                count=count+1;
            end
        end
    end
            t(z)=count;
            count=0;
end


t1=1:256;
n1=0:255;
count=0;
for z=1:256
    for i=1:M
        for j=1:N           
            if Green(i,j)==z-1
                count=count+1;
            end
        end
    end
            t1(z)=count;
            count=0;
end

t2=1:256;
n=0:255;
count=0;
for z=1:256
    for i=1:M
        for j=1:N        
            if Blue(i,j)==z-1
                count=count+1;
            end
        end
    end
            t2(z)=count;
            count=0;
end


axes(handles.axes2);
plot(n,t,'Red', n, t1, 'Green', n, t2, 'Blue');
xlabel('Pixel Intensity');
ylabel('Number of Pixels');
xlim([0 256])
end



function angle_rotate_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function angle_rotate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function scale_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function scale_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scale2_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function scale2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rescalebutton.
function rescalebutton_Callback(hObject, eventdata, handles)

global im2
inputImage=im2;

% scale = input('scale = ')  ;           %# The resolution scale factors: [rows columns]
scale(1)=str2num(get(handles.scale,'String'));
scale(2)=str2num(get(handles.scale2,'String'));
oldSize = size(inputImage);                   %# Get the size of your image
newSize = max(floor(scale.*oldSize(1:2)),1);  %# Compute the new image size

%# Compute an upsampled set of indices:

rowIndex = min(round(((1:newSize(1))-0.5)./scale(1)+0.5),oldSize(1));
colIndex = min(round(((1:newSize(2))-0.5)./scale(2)+0.5),oldSize(2));

%# Index old image to get new image:

outputImage = inputImage(rowIndex,colIndex,:);


%axes(handles.axes2);
figure()
imshow(inputImage);
figure()  
imshow(outputImage);


% --- Executes on button press in anglebutton_90.
function anglebutton_90_Callback(hObject, eventdata, handles)

global im2
image=im2; 

angle = 90;
outputImage = rotating(image,angle);

axes(handles.axes2)
imshow(outputImage)


% --- Executes on button press in anglebutton_180.
function anglebutton_180_Callback(hObject, eventdata, handles)

global im2
image=im2; 

angle=180;
outputImage = rotating(image,angle);

axes(handles.axes2)
imshow(outputImage)


% --- Executes on button press in anglebutton_270.
function anglebutton_270_Callback(hObject, eventdata, handles)

global im2
image=im2; 

angle=270;
outputImage = rotating(image,angle);

axes(handles.axes2)
imshow(outputImage)


% --- Executes on button press in anglebutton.
function anglebutton_Callback(hObject, eventdata, handles)

global im2
image=im2; 

angle=str2double(get(handles.angle_rotate,'String'));
outputImage = rotating(image,angle);

axes(handles.axes2)
imshow(outputImage)


% --- Executes on selection change in mask_size.
function mask_size_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function mask_size_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sigma_value_Callback(hObject, eventdata, handles)

sigma = (get(handles.sigma_value,'Value'));
set(handles.sigma_value_display,'String',num2str(sigma,2))


% --- Executes during object creation, after setting all properties.
function sigma_value_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function sigma_value_display_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function sigma_value_display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in method_selection.
function method_selection_Callback(hObject, eventdata, handles)

method_number = get(handles.method_selection,'Value');


%For Gaussian Blur 
if method_number == 2;
    %check_string = get(handles.text2,'String');
   % if strcmp(check_string,'D Value')==1
    set(handles.sigma_text,'String','Sigma  ')
    set(handles.masksize_text,'String','Mask Size   ')
    set(handles.sigma_value,'Max',5,'Min',0.1,'Value',0.1,'String','0.1')
    set(handles.sigma_value_display,'Max',5,'Min',0.1,'Value',0.1,'String','0.1')
    set(handles.mask_size,'Style','popupmenu','String',['2x2'; '3x3'; '4x4'; '5x5'])
   % end
end



% --- Executes during object creation, after setting all properties.
function method_selection_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in blurringbutton.
function blurringbutton_Callback(hObject, eventdata, handles)


global im2
im_original = im2;

method_number = get(handles.method_selection,'Value'); % selected method
if method_number == 2
    mask_size = get(handles.mask_size,'String');       % size of the mask
    mask_size_number = get(handles.mask_size,'Value'); % size of the mask
    if mask_size_number == 1;
        n = 2;
        m = n;
    elseif mask_size_number  == 2;
        n = 3;
        m = n;
    elseif mask_size_number  == 3;
        n = 4;
        m = n;
    elseif mask_size_number  == 4;
        n = 5;
        m = n;
    end
    sigma = get(handles.sigma_value,'Value');                     % read value of the sigma
    set(handles.sigma_value_display,'String',num2str(sigma,2));     % set value of sigma for diaplay
    if method_number ==2,
        im_filtered = image_blur(im_original,sigma,n,m);          % filter image using image_blur function
    end
end
    
    axes(handles.axes2);
    imshow(im_filtered);
    
     % Function to calculate filtered gaussian blur image
function res = image_blur(data,sigma,n,m)

x = double(data);                       % convert from uint8 to double
filter1 = i_gauss(n,sigma);             % applying i_gauss function to find filtered image
% % center = (n+1)/2;
% % res = zeros(size(data));
% % for k=1:n
% %    ki = k-center
% %   for l=1:n
% %       li = l-center
% %       res(center:end-center+1,center:end-center+1) = ...
% %           res(center:end-center+1,center:end-center+1) + ...
% %           x(center-ki:end-ki-center+1,center-li:end-center-li+1).*filter1(k,l);
% %   end
% % end 
% % res = uint8(res);

res = uint8(convn(x,filter1,'same'));   % convolve and convert to uint8 from double 

%----------------------------------------------------------
% Function to calculate gaussian mask matrix
function res = i_gauss(x_dim,sigma)
% Sigma is the standard deviation of the gaussian function
n = x_dim;           % Dimensions of gaussian mask matrix i.e., nxn
for i = 1:n,
    for j = 1:n
        M_temp = [j-(n+1)/2 i-(n+1)/2]';       
        M(i,j) = gauss(M_temp(1),sigma)*gauss(M_temp(2),sigma); % mask
    end
end
res = M/sum(sum(M)); % to normalize the mask

%----------------------------------------------------------
%Function to calculate Gaussian function
function res = gauss(x,sigma)

res = exp(-x^2/(2*sigma^2))/(sigma*sqrt(2*pi));
    

%----------------------------------------------------------




function threshold_value_display_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function threshold_value_display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function threshold_value_Callback(hObject, eventdata, handles)


threshold = (get(handles.threshold_value,'Value'));
set(handles.threshold_value_display,'String',threshold)

% --- Executes during object creation, after setting all properties.
function threshold_value_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in thresholdbutton.
function thresholdbutton_Callback(hObject, eventdata, handles)

global im2;
image = im2;
image_black = toGrayscale(image);
[i, j] = size(image_black);

threshold = (get(handles.threshold_value,'Value'));
THR=threshold;

for M=1:i
    for N=1:j
        if image_black(M,N) >= THR;
            image_black(M,N)=255;
        else
            image_black(M,N)=0;
        end
    end
end


axes(handles.axes2);
imshow(image_black)



% --- Executes on button press in egdedetectionbutton.
function egdedetectionbutton_Callback(hObject, eventdata, handles)


global im2

I_gray=double(toGrayscale(im2));
[w,h]=size(I_gray);

x_kernel=[-1,0,1;-2,0,+2;-1,0,+1];
y_kernel=[-1,-2,-1;0,0,0;1,2,1];

I_new_x=zeros(w,h);
I_new_y=zeros(w,h);

for i=2:w-1
    for j=2:h-1
        I_new_x(i,j)=I_gray(i,j)*x_kernel(2,2)+I_gray(i-1,j-1)*x_kernel(1,1)+I_gray(i+1,j+1)*x_kernel(3,3)...
          +I_gray(i,j+1)*x_kernel(2,3)+ I_gray(i-1,j)*x_kernel(1,2)+I_gray(i,j-1)*x_kernel(2,1)...
          +I_gray(i-1,j+1)*x_kernel(1,3)+I_gray(i+1,j-1)*x_kernel(3,1)+I_gray(i+1,j)*x_kernel(3,2);
      
        I_new_y(i,j)=I_gray(i,j)*y_kernel(2,2)+I_gray(i-1,j-1)*y_kernel(1,1)+I_gray(i+1,j+1)*y_kernel(3,3)...
            +I_gray(i,j+1)*y_kernel(2,3)+ I_gray(i-1,j)*y_kernel(1,2)+I_gray(i,j-1)*y_kernel(2,1)...
            +I_gray(i-1,j+1)*y_kernel(1,3)+I_gray(i+1,j-1)*y_kernel(3,1)+I_gray(i+1,j)*y_kernel(3,2);
    end
end
magnitude=sqrt(I_new_y.^2+I_new_x.^2);

edgeImage = magnitude > 240;
axes(handles.axes2);
imshow(edgeImage); 


% --- Executes on button press in grayscalebutton.
function grayscalebutton_Callback(hObject, eventdata, handles)

global im2
returnedImage=toGrayscale(im2);
axes(handles.axes2);
imshow(returnedImage);



% --------------------------------------------------------------------
function Edit_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Reset_Callback(hObject, eventdata, handles)

axes(handles.axes2);
hold off;
cla reset;

set(handles.static,'String','');
set(handles.axes1,'xtick',[],'ytick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
set(handles.axes2,'xtick',[],'ytick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);

% --------------------------------------------------------------------
function Clear_Callback(hObject, eventdata, handles)

axes(handles.axes1);
hold off;
cla reset;
axes(handles.axes2);
hold off;
cla reset;
set(handles.static,'String','');
set(handles.axes1,'xtick',[],'ytick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
set(handles.axes2,'xtick',[],'ytick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);



%----------------------------------------------------------------------
function outputImage = rotating(image,angle)

[rowsi,colsi,z]= size(image); 


rads=-(2*pi*angle/360);  

%calculating array dimesions such that  rotated image gets fit in it exactly.
% we are using absolute so that we get  positve value in any case ie.,any quadrant.

rowsf=ceil(rowsi*abs(cos(rads))+colsi*abs(sin(rads)));                      
colsf=ceil(rowsi*abs(sin(rads))+colsi*abs(cos(rads)));                     

% define an array withcalculated dimensionsand fill the array  with zeros ie.,black
outputImage=uint8(zeros([rowsf colsf 3 ]));

%calculating center of original and final image
xo=ceil(rowsi/2);                                                            
yo=ceil(colsi/2);

midx=ceil((size(outputImage,1))/2);
midy=ceil((size(outputImage,2))/2);

% in this loop we calculate corresponding coordinates of pixel of A 
% for each pixel of C, and its intensity will be  assigned after checking
% weather it lie in the bound of A (original image)
    for i=1:size(outputImage,1)
        for j=1:size(outputImage,2)                                                       

             x= (i-midx)*cos(rads)+(j-midy)*sin(rads);                                       
             y= -(i-midx)*sin(rads)+(j-midy)*cos(rads);                             
             x=round(x)+xo;
             y=round(y)+yo;

             if (x>=1 && y>=1 && x<=size(image,1) &&  y<=size(image,2) ) 
                  outputImage(i,j,:)=image(x,y,:);  
             end
         end
    end
