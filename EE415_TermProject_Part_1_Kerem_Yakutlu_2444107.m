%EE415 - Term Project - Part 1
%Kerem Yakutlu - 2444107
%% 
% *Test Image Projecfion*
% 
% In this code block, I have written how to project the test image.

I = cell2mat(struct2cell(load('square.mat'))); %Loading the image, as an input.

Beams = 100; %Number of beams as an input.
Step_Size = deg2rad(1); %Specifying beams and step size, as an input.

Projection_of_the_Test_Image = Projector(I, Beams, Step_Size); %Using the projecting function for the test image.

%%%Plotting.

Plotting_the_Projections(Projection_of_the_Test_Image, Beams, Step_Size);
%% 
% *Complex Image Projecfion*
% 
% In this code block, I have written how to project a more complex image.


I = cell2mat(struct2cell(load('lena.mat'))); %Loading the image, as an input.

Beams = 100; %Number of beams as an input.
Step_Size = deg2rad(1); %Specifying beams and step size, as an input.

Projection_of_the_Complex_Image = Projector(I, Beams, Step_Size); %Using the projecting function for the test image.

%%%Plotting.

Plotting_the_Projections(Projection_of_the_Complex_Image, Beams, Step_Size);
%% 
% *Function for the Projection*
% 
% These are the functions for projection and plotting the sample images.

function [Projection_of_the_Image] = Projector(I, Beams, Step_Size) %Getting the image by taking the .mat form.
    [M , N] = size(I); %Taking the size of I, where M is the number of rows and N is the number of columns.
    
    t = linspace(-M/sqrt(2),M/sqrt(2),Beams); %Specifying t and theta by using linspace.
    theta = (0:Step_Size:(pi - Step_Size));
    
    Result_Matrix = zeros(length(theta), length(t)); %Obtaining the result matrix, which will be used after the back projection matrix.
    
    X_Index_Vector = (-M/2):(M/2); %Setting x vector for indexing.
    Y_Index_Vector = (-M/2):(M/2); %Setting y vector for indexing.
    
    %%%Calculating the projection and backprojection integrals.
    
    for i = 1:length(t)
        for j = 1:length(theta)  
            Points_of_X = X_Index_Vector; %Straightforward appointing X values.
            Points_of_Y = ( t(i) - X_Index_Vector*cos(theta(j)) ) / sin(theta(j)); %Calculating Y intersections for integer Y values
    
            Points_of_X(Points_of_Y < -M/2) = []; %Clearing the low values, lower than M/2
            Points_of_Y(Points_of_Y < -M/2) = [];
    
            Points_of_X(Points_of_Y > M/2) = []; %Clearing the high values, higher than M/2
            Points_of_Y(Points_of_Y > M/2) = [];
    
            Combined_Vectors = [Points_of_X' Points_of_Y']; %Combining the vectors of X and Y.
    
            %%%%%%%%%%%%%%%%%%%%%%%%
    
            Points_of_Y = Y_Index_Vector; %Straightforward appointing X values.
            Points_of_X = ( t(i) - Y_Index_Vector*sin(theta(j)) ) / cos(theta(j)); %Calculating X intersections for integer X values
    
            Points_of_Y(Points_of_X < -M/2) = []; %Clearing the low values, lower than M/2
            Points_of_X(Points_of_X < -M/2) = [];
    
            Points_of_Y(Points_of_X > M/2) = []; %Clearing the high values, higher than M/2
            Points_of_X(Points_of_X > M/2) = [];
    
            Recombined_Vector_for_Points_of_X = [Combined_Vectors(:,1);Points_of_X']; %Combining another set of vectors.
            Recombined_Vector_for_Points_of_Y = [Combined_Vectors(:,2);Points_of_Y'];
    
            Combined_Vectors = [Recombined_Vector_for_Points_of_X Recombined_Vector_for_Points_of_Y]; %Concentration of the vectors.
    
            %%%%%
    
            Coordinates_Sorted = unique(round(sortrows(Combined_Vectors,1),4),'rows'); %Row sorting for consecutive intersection points.
            
            Distance = zeros(1,length(Coordinates_Sorted)-1); %Creating distance vector, for calculation of projection (Projection Result = Square Value * Distance between Intersection Points) 
            Midpoints_X = zeros(1,length(Coordinates_Sorted)-1); %Creating midpoint X and Y vectors
            Midpoints_Y = zeros(1,length(Coordinates_Sorted)-1);
    
            for k = 1:(length(Coordinates_Sorted(:,1)) -1)
    
                if (length(Coordinates_Sorted) == 1) %If there is only one coordinate set, the distance would be nothing but 0.
                    Distance(k) = 0;
                    Midpoints_X(k) = 1;
                    Midpoints_Y(k) = 1;
                else
                    Distance(k) = sqrt((Coordinates_Sorted(k,1)-Coordinates_Sorted(k+1,1))^2+(Coordinates_Sorted(k,2)-Coordinates_Sorted(k+1,2))^2); %Finding the distance between consecutive intersection points.
                    Midpoints_X(k) = (Coordinates_Sorted(k+1,1) + Coordinates_Sorted(k,1))/2;
                    Midpoints_Y(k) = (Coordinates_Sorted(k+1,2) + Coordinates_Sorted(k,2))/2; %Finding the midpoints
                end
            end
    
            Row_Data = M/2 - floor(Midpoints_Y);
            Column_Data = M/2 + ceil(Midpoints_X); %Row and column data must be positive for the sake of MATLAB indices. 
    
                Distance(isinf(Distance)) = []; %Clearing the infinite elements, for the easeness of calculations.
                Row_Data(isinf(Row_Data)) = [];
                Column_Data(isinf(Column_Data)) = [];
             
                Distance(isnan(Distance)) = []; %Clearing the NaN elements, for the easeness of calculations.
                Row_Data(isnan(Row_Data)) = [];
                Column_Data(isnan(Column_Data)) = [];
    
                for n = 1:length(Row_Data)
                    Result_Matrix(j,i) = Result_Matrix(j,i) + Distance(n) * I(Row_Data(n), Column_Data(n)); %Taking the line integral of projections. (Projection Result = Square Value * Distance between Intersection Points) 
                end
    
        end

        Projection_of_the_Image = Result_Matrix;
    
end
end

function Plotting_the_Projections(Projection_of_The_Image, Beams, Step_Size) %Function for plotting the image.
    Place_Counter = 1; %Counter will be used in the for loop.
    figure;
    hold on;
    
    Theta_for_Plotting = (0:Step_Size:(pi - Step_Size)); %Theta angles are determined for the plotting.

    sgtitle("Projection Plots"); %Main title for the subplots.
    
    for m = round(1:length(Theta_for_Plotting)/12:length(Theta_for_Plotting))
        subplot(3,4,Place_Counter); %Subplots are used.
        plot(1:Beams,Projection_of_The_Image(m,:)) %Plotting every projection.
        title(num2str(rad2deg(Theta_for_Plotting(m))) + "°", FontSize=8); %Implementing the title.
        xlabel('Beam Numbers', FontSize=6); %Implementing the x label.
        ylabel('Projection Data', FontSize=6); %Implementing the y label.
        Place_Counter = Place_Counter + 1; %Counter is increased by 1 at every for loop iteration.
    end
    hold off;
end