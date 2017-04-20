clear all;
close all;
clc;

%% Fast Fratcal Coding Method for the detection of microcalcificationin mammograms
% Author 1: Mayank Arya, admin@mamexo.me
% Author 2: Ram Lakhan , lakhan.iitp@gmail.com
% Implementation in brief
% =======================
% In this implementaiton of microcalcification detection, for a mammogram
% medical image range and domain blocks are classified and two adaptive
% thresholds T1 and T2 are are used for the binaryzation before which the
% restored image is median filtered.
% Input
% ======
% A mammogram image name which is need to be set in the variable @imagename
% Output
% ======
% The input mammogram and the dectected mamaogram are shown in two figures
% How to use
% ===========
% To use this program follow the following steps
% STEP 1: Set the image name same as the name of the input image
% STEP 2: Run the program
% STEP 3: Output windows will popup automatically after processing is
% complete
% The following program has been tested on these files from the mias
% database
% imagenames: 'mdb008.pgm','mdb028.pgm','mdb134.pgm','mdb184.pgm'
% Variables used 
% ===============
% 

%% STEP 1 
% Set the input image name here
imagename = 'mdb184.pgm';


%% image is read and resized to size of interest and image paritioning size i.e size of range block and domain block is defined here.
image = imread(imagename);
resized_image = imresize(image,[32,32]);
[l, b] = size(resized_image);
show_image = resized_image;
input_image = double(resized_image);
modeled_input_image = double(zeros(l));
rangesize = 4;
domainsize = 8;

%% Encoding And Decoding of input image
% Search for a matching domain block D which has similar structures as that of R.
% Each of the blocks in the domain search pool is scaled to the size of range blocks to encode range blocks. It is then
% compared with i R using isometric transformations and intensity offset as well as intensity scaling (contrast scaling).
% The transformation w for the range blocks R is defined by:
%                           w(R) = s ×(D) + o
% where s and o are intensity scaling and offset factor.
% s and o are given as 
% s = (N*sum(sum(R.D))-(sum(sum(R)) * sum(sum(D))))/(N*sum(sum(R^.2)) - sum(sum(D.^2);     
% o = (sum(sum(D)) - s_i * sum(sum(R)))/N;
% where N is total number of pixels in the range block and D is the spatial
% contracted domain blocks under following isometric transformations
% (1) identity, 
% (2) rotation through +90°, 
% (3) rotation through +180°, 
% (4) rotation through ?90°, 
% (5) reflection about mid-vertical axis, 
% (6) reflection about mid-horizontal axis, 
% (7) reflection about first diagonal, and 
% (8) reflection about second diagonal.
% The domain block which minimizes the distortion i e is chosen. The distortion is expressed as:
%  error = sum(sum((R - (s_i * D+o_i)).^2))





for rangeblock_x = 1:rangesize:b-rangesize+1
    for rangeblock_y = 1:rangesize:b-rangesize+1
       
        %% Calculation of range block
        range_block = input_image(rangeblock_x:rangeblock_x + rangesize-1,rangeblock_y:rangeblock_y + rangesize-1); %range block matrix
        
        
        %% Calculation of Dynamic range for range blocks
        min_value = min(min(range_block));
        max_value = max(max(range_block));
        modeled_input_image(rangeblock_x:rangeblock_x + rangesize-1,rangeblock_y:rangeblock_y + rangesize-1) = range_block;
        dynamic_range_range = 1- double(min_value)/double(max_value);   
        
        %% Test for dynamic range of range blocks to be between 0.05 to 1 for encoded by fractal encoding method
        if dynamic_range_range > 0.05 && dynamic_range_range <= 1 
            
            min_error = Inf;
            
            for domainblock_x= 1:4:b-domainsize+1
                for domainblock_y= 1:4:b-domainsize+1

                    %% Calculation of domain blocks
                    domain_block = input_image(domainblock_x:domainblock_x + domainsize-1,domainblock_y:domainblock_y + domainsize-1);
                    
                    %% Calcultaion of dynamic range for domain blocks
                    min_value = min(min(domain_block));
                    max_value = max(max(domain_block));
                    dynamic_range_domain = 1- double(min_value)/double(max_value);
                    
                    %% Test for dynamic range of domian blocks for to be used in the encoding of range blocks
                    if dynamic_range_domain >= 0.07 && dynamic_range_domain <= 0.35   
                      
                        %% Calculation of 8 affine transforamtion to be applied on the domain blocks for finding the best mapping between domain and range blocks
                        % Resizing the domain block to the size of the
                        % range blocks
                        D_scale = imresize(domain_block,[rangesize,rangesize]);
                        
                        Domain = D_scale;
                        % Transformation 1: Identity 
                        DT_MAT(1:rangesize, 1:rangesize, 1) = Domain;
                        % Transformation 2: Rotaion by +90
                        DT_MAT(1:rangesize, 1:rangesize, 2) = rot90(Domain);
                        % Transformation 3: Rotation by +180
                        DT_MAT(1:rangesize, 1:rangesize, 3) = rot90(rot90(Domain));
                        % Transformation 4: Rotation by -90
                        DT_MAT(1:rangesize, 1:rangesize, 4) = rot90(rot90(rot90(Domain)));
                        % Transformation 5: Reflection about mid horizontal
                        % axis
                        DT_MAT(1:rangesize, 1:rangesize, 5) = flipud(Domain);
                        % Transformation 6: Reflection about mid vertical
                        % axis
                        DT_MAT(1:rangesize, 1:rangesize, 6) = fliplr(Domain);
                        % Transformation 7: Rotation about first diagonal
                        DT_MAT(1:rangesize, 1:rangesize, 7) = transpose(Domain);
                        % Transformation 8: Rotation about second diagonal
                        DT_MAT(1:rangesize, 1:rangesize, 8) = rot90(rot90(transpose(Domain)));
                        % N is the number of pixels in the range block
                        N = rangesize * rangesize;
                        
                        %% Finding the best mapping for the range block using current domain block and its 7 other affine transformation
                        for Loop=1:8
                            Domain = DT_MAT(1:rangesize, 1:rangesize, Loop);
                            %% Calulation for intensity scaling factor and intensity offset factor
                            s_i = (N*sum(sum(range_block.*Domain))-(sum(sum(range_block)) * sum(sum(Domain))))/(N*sum(sum(range_block.^2)) - sum(sum(range_block))^2);     
                            o_i = (sum(sum(Domain)) - s_i * sum(sum(range_block)))/N;
                            %% Calculation of the error between the range block and the domain block's transformation if this error is less than the minimum error then use this in the decoding process.
                            error = sum(sum((range_block - (s_i * Domain+o_i)).^2));
                            if error < min_error
                                min_error = error;
                                modeled_input_image(rangeblock_x:rangeblock_x + rangesize-1,rangeblock_y:rangeblock_y + rangesize-1) = s_i *Domain + o_i;
                            end
                        end
                    end
               end
           end
        end
    end
end
%% Finding the microcalcification by finding the difference between the modelled image and the original image.

% Renaming the variables for better readibilty
original_image = double(input_image);
restored_image = double(modeled_input_image);


%% Applying the median filter to the modelled image or restored image 
med_fil_res_image = medfilt2(restored_image);

%% Finding the difference between the original image and the median filtered modelled image 
diff_image = double(original_image) - double(med_fil_res_image);

%% Calculation of two adpative thresholds T1 and T2

% finding the histogram of the gray level difference image 
hist_count = histc(diff_image(:),0:255);
hist_data = zeros(size(hist_count),2);
hist_data(:,1) = [0:255];
hist_data(:,2) = [hist_count+1];


%% Uncomment the following code if you want to see the gray level histogram and the see the calculated centroids for the gray level difference image
% figure;
% plot(hist_data(:,1),hist_data(:,2),'.');
% axis([0 255 0 250]);
% title 'Histogram of difference image';

%CALCULATION pf histogram centres using k means algorithm 
opts = statset('Display','final');
[idx,C] = kmeans(hist_data,2,'Distance','cityblock','Replicates',5,'Options',opts);

% figure;
% plot(hist_data(idx==1,1),hist_data(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(hist_data(idx==2,1),hist_data(idx==2,2),'b.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Centroids','Location','NW')
% title ('Cluster Assignments and Centroids')
% hold off


%% Calculation of T1 and T2 using the first cluster center and the mean and standard devaition of the image respectively.

%% Calculation of T1
C = sort(C);
t = 0.5; % t is constant
T1 = t * C(1);

%% Calcutaion of T2
m = mean(mean(original_image));
sigma = std(std(original_image));
k = 1.56;
T2 = m - k * sigma;


%% Binarization process

%create a variable to hold the image after binarization process
final_image = zeros(size(restored_image));

for input_image=1:l
    for y=1:l
        % Test if the pixel value in the difference image is greater than threshold T1 and the gray level value in the original image is greater than the threshold T2. 
        if diff_image(input_image,y) > T1 && original_image(input_image,y) > T2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
            final_image(input_image,y) = 255;
        else
            final_image(input_image,y) = 0;
        end
    end
end

%% Showing the input image or mammogram and the processed image afer binarization or the detected microcalcification
figure;
imshow([show_image final_image]);

%% Uncoment the upcoming section to see Some intermediate results.
% figure;
% imshow(show_image);
% title('Input mammogram image')
% 
% figure;
% imshow(restored_image);
% title('Modeled input mammogram using fratal modelling')
% 
% figure;
% imshow(med_fil_res_image);
% title('Median filtered modelled image')
% 
% 
% figure;
% imshow(diff_image);
% title('Detected microcalcification by difference mechanism between original and modelled image')
% 
% figure;
% imshow(final_image);
% title('Detected microcalcification after applying bianrization process ')