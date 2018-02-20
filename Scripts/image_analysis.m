%%% This script uses the map and digital elevation model of the Domboshawa region
%%% to incorporate elevation information in the distribution of the present lichen species.
%%% No statistical analysis is included in the script, but was conducted with R.
%%% The "%%" symbols indicate the beginning of a new section.


cd 'working_directory' % set working directory
cd
%% create greyscale elevation heatmap

resolution = 256;
cmap = jet(resolution);
x = cmap(:,1);
y = cmap(:,2);
z = cmap(:,3);
intensity = gray(256);
intensity = intensity(:,1);
fx = scatteredInterpolant(x,y,z,intensity,'nearest'); % correlates jet map values with greyscale values
img = double(imread('data_path'))/255; % import digital elevation file
imgOut = fx(img(:,:,1),img(:,:,2),img(:,:,3)); % applies greyscale map to image
imshow(imgOut)
%% add elevation values to RGB image

RGB= imread('data_path'); % import map file � slightly higher resolution than heatmap image
SizeRGB= size(RGB)
RGB= RGB(:,:,1:3);
A= im2double(RGB); % create copy in double format to add elevation matrix later (4th dimension is binary in uint8 images)
imshow(RGB)
Heat_grey=double(imgOut); % values range from 0 to 1 (low to high elevation)
Heat_grey= imresize(Heat_grey, [SizeRGB(1), SizeRGB(2)], 'nearest'); % resize heat through nearest neighbor interpolation (multiplies pixels)
% some elevation information will be artificial, but only a very small part
size(Heat_grey)==size(RGB(:,:,1)) % chech if rgb-heatmap matrices have the same size
RGB_map= cat(3, A(:,:,1), A(:,:,2), A(:,:,3), Heat_grey); % concatenate matrices on third dimension
size(RGB_map) % RGB_map stores RGB values and elevation information (the elevation matrix is mostly important to append to the cluster matrix)
%% crop image 
% produce equal size images, which include all species and show altitudinal gradient
imtool(RGB)
R1= RGB(2369:3343, 548:1758, 1:3); % isolate first area
imshow(R1)
R2= RGB(1362:2336 , 260:1470, 1:3);
imshow(R2)
R3= RGB(20:994, 372:1582, 1:3);
imshow(R3)
R4= RGB(43:1017, 2501:3711, 1:3);
imshow(R4)
R5= RGB(1362:2336, 1842:3052, 1:3);
imshow(R5)
R6= RGB(1146:2120, 3468:4678, 1:3);
imshow(R6)
% size= 975x1211x3
%% find validation points per species
imtool(R1)
R1_red= R1(331:337, 583:598, 1:3)
R1_pink= R1(491:504, 207:226, 1:3)
R1_green= R1(472:477, 1095:1112, 1:3)
R1_black= R1(502:512, 1031:1046, 1:3)
R1_stone= R1(317:328, 508:519, 1:3)

imtool(R2)
R2_red= R2(210:216, 301:310, 1:3)
R2_pink= R2(512:524, 938:949, 1:3)
R2_green= R2(720:725, 653:663, 1:3)
R2_black= R2(536:549, 202:219, 1:3)

imtool(R3)
R3_red=  R3(516:524, 270:285, 1:3)
R3_pink= R3(584:594, 1066:1077, 1:3)
R3_green= R3(416:424, 675:685, 1:3)
R3_black= R3(628:639, 636:645, 1:3)

imtool(R4)
R4_red= R4(461:470, 534:548, 1:3)
R4_pink= R4(580:591, 835:848, 1:3)
R4_green= R4(309:316, 814:823, 1:3)
R4_black= R4(252:262, 353:365, 1:3)

imtool(R5)
R5_red= R5(167:175, 378:391, 1:3)
R5_pink= R5(522:531, 558:567, 1:3)
R5_green= R5(676:686, 1046:1053, 1:3)
R5_black= R5(675:682, 333:348, 1:3)

imtool(R6)
R6_red= R6(219:226, 806:813, 1:3)
R6_pink= R6(508:515, 226:239, 1:3)
R6_green= R6(344:351, 1112:1121, 1:3)
R6_black= R6(500:507, 492:504, 1:3)

%% k-means clustering
% cluster indices may correspond to different colours in each run, so values below may have to change every time 
% check which colour corresponds to each index each time and use the appropriate values
Size1= size(R1);
I_data= double(R1); % convert to double class
I_data= reshape(I_data, [Size1(1)*Size1(2) 3]); 
[idx centroids]= kmeans(I_data, 3, 'Replicates', 10);  
I_cluster1= reshape(idx, [Size1(1) Size1(2)]); 
imtool(I_cluster1, []) 
Group1= find(I_cluster1==1); % number of pixels belonging to cluster 1 � blacl
Proportion_1= size(Group1)/numel(I_cluster1);
Validation_1= size(find(I_cluster1(502:512, 1031:1046)==1))/numel(I_cluster1(502:512, 1031:1046))
Group2= find(I_cluster1==2); % pink
Proportion_2= size(Group2)/numel(I_cluster1); % 27.93% 
Validation_2= size(find(I_cluster1(491:504, 207:226)==2))/numel(I_cluster1(491:504, 207:226))
Group3= find(I_cluster1==3); % red
Proportion_3= size(Group3)/numel(I_cluster1); % 26.59
Validation_3= size(find(I_cluster1(331:337, 583:598)==3))/numel(I_cluster1(331:337, 583:598)) 
Proportion_1(1)+Proportion_2(1)+Proportion_3(1) % sum equals 1
Validation1_mean= (Validation_1(1)+Validation_2(1)+Validation_3(1))/3 % 98.99%

%R2
Size2= size(R2);
I_data= double(R2); % convert to double class
I_data= reshape(I_data, [Size2(1)*Size2(2) 3]); 
[idx centroids]= kmeans(I_data, 3, 'Replicates', 10);  
I_cluster2= reshape(idx, [Size2(1) Size2(2)]); 
imtool(I_cluster2, []) 
Group1= find(I_cluster2==1); % number of pixels belonging to cluster 1 � pink (40.53%)
Proportion_1= size(Group1)/numel(I_cluster2); % proportion of pixels belonging to cluster 1
Validation_1= size(find(I_cluster2(512:524, 938:949)==1))/numel(I_cluster2(512:524, 938:949))
Group2= find(I_cluster2==2); % red
Proportion_2= size(Group2)/numel(I_cluster2); % 44.38%
Validation_2= size(find(I_cluster2(210:216, 301:310)==2))/numel(I_cluster2(210:216, 301:310))
Group3= find(I_cluster2==3); % black
Proportion_3= size(Group3)/numel(I_cluster2); % 15.10
Validation_3= size(find(I_cluster2(536:549, 202:219)==3))/numel(I_cluster2(536:549, 202:219)) 
Proportion_1(1)+Proportion_2(1)+Proportion_3(1) % sum equals 1
Validation1_mean= (Validation_1(1)+Validation_2(1)+Validation_3(1))/3 % 92.38%

%R3
Size3= size(R3);
I_data= double(R3); % convert to double class
I_data= reshape(I_data, [Size2(1)*Size2(2) 3]); 
[idx centroids]= kmeans(I_data, 3, 'Replicates', 10);  
I_cluster3= reshape(idx, [Size2(1) Size2(2)]); 
imtool(I_cluster3, []) 
Group1= find(I_cluster3==1); % number of pixels belonging to cluster 1 � red 
Proportion_1= size(Group1)/numel(I_cluster3); % proportion of pixels belonging to cluster 1 (43.62%)
Validation_1= size(find(I_cluster3(516:524, 270:285)==1))/numel(I_cluster3(516:524, 270:285))
Group2= find(I_cluster3==2); % black
Proportion_2= size(Group2)/numel(I_cluster3); % 14.4%
Validation_2= size(find(I_cluster3(628:639, 636:645)==2))/numel(I_cluster3(628:639, 636:645))
Group3= find(I_cluster3==3); % pink
Proportion_3= size(Group3)/numel(I_cluster3); % 41.98%
Validation_3= size(find(I_cluster3(580:591, 835:848)==3))/numel(I_cluster3(580:591, 835:848)) 
Proportion_1(1)+Proportion_2(1)+Proportion_3(1) % sum equals 1
Validation1_mean= (Validation_1(1)+Validation_2(1)+Validation_3(1))/3 % 57.18%

%R4
Size4= size(R4);
I_data= double(R4); % convert to double class
I_data= reshape(I_data, [Size2(1)*Size2(2) 3]); 
[idx centroids]= kmeans(I_data, 3, 'Replicates', 10);  
I_cluster4= reshape(idx, [Size2(1) Size2(2)]); 
imtool(I_cluster4, []) 
Group1= find(I_cluster4==1); % number of pixels belonging to cluster 1 � pink
Proportion_1= size(Group1)/numel(I_cluster4); % proportion of pixels belonging to cluster 1 (42.31%)
Validation_1= size(find(I_cluster4(580:591, 835:848)==1))/numel(I_cluster4(580:591, 835:848))
Group2= find(I_cluster4==2); % black
Proportion_2= size(Group2)/numel(I_cluster4); %13.53
Validation_2= size(find(I_cluster4(252:262, 353:365)==2))/numel(I_cluster4(252:262, 353:365))
Group3= find(I_cluster4==3); % red
Proportion_3= size(Group3)/numel(I_cluster4); % 44.17
Validation_3= size(find(I_cluster4(461:470, 534:548)==3))/numel(I_cluster4(461:470, 534:548)) 
Proportion_1(1)+Proportion_2(1)+Proportion_3(1) % sum equals 1
Validation1_mean= (Validation_1(1)+Validation_2(1)+Validation_3(1))/3 % 84.67%

%R5
Size5= size(R5);
I_data= double(R5); % convert to double class
I_data= reshape(I_data, [Size2(1)*Size2(2) 3]); 
[idx centroids]= kmeans(I_data, 3, 'Replicates', 10);  
I_cluster5= reshape(idx, [Size2(1) Size2(2)]); 
imtool(I_cluster5, []) 
Group1= find(I_cluster5==1); % number of pixels belonging to cluster 1 � pink
Proportion_1= size(Group1)/numel(I_cluster5); % proportion of pixels belonging to cluster 1 (44.25%)
Validation_1= size(find(I_cluster5(522:531, 558:567)==1))/numel(I_cluster5(522:531, 558:567))
Group2= find(I_cluster5==2); % black
Proportion_2= size(Group2)/numel(I_cluster5); % 12.71%
Validation_2= size(find(I_cluster5(675:682, 333:348)==2))/numel(I_cluster5(675:682, 333:348))
Group3= find(I_cluster5==3); % red
Proportion_3= size(Group3)/numel(I_cluster5); % 43.05%
Validation_3= size(find(I_cluster5(461:470, 534:548)==3))/numel(I_cluster5(461:470, 534:548)) 
Proportion_1(1)+Proportion_2(1)+Proportion_3(1) % sum equals 1
Validation1_mean= (Validation_1(1)+Validation_2(1)+Validation_3(1))/3 % 76.44%

%R6
Size6= size(R6);
I_data= double(R6); % convert to double class
I_data= reshape(I_data, [Size2(1)*Size2(2) 3]); 
[idx centroids]= kmeans(I_data, 3, 'Replicates', 10);  
I_cluster6= reshape(idx, [Size2(1) Size2(2)]); 
imtool(I_cluster6, []) 
Group1= find(I_cluster6==1); % number of pixels belonging to cluster 1 � black
Proportion_1= size(Group1)/numel(I_cluster6); % proportion of pixels belonging to cluster 1 (11.98%)
Validation_1= size(find(I_cluster6(500:507, 492:504)==1))/numel(I_cluster6(500:507, 492:504))
Group2= find(I_cluster6==2); % pink
Proportion_2= size(Group2)/numel(I_cluster6); % 44.45%
Validation_2= size(find(I_cluster6(508:515, 226:239)==2))/numel(I_cluster6(508:515, 226:239))
Group3= find(I_cluster6==3); % red
Proportion_3= size(Group3)/numel(I_cluster6); % 43.57%
Validation_3= size(find(I_cluster6(219:226, 806:813)==3))/numel(I_cluster6(219:226, 806:813)) 
Proportion_1(1)+Proportion_2(1)+Proportion_3(1) % sum equals 1
Validation1_mean= (Validation_1(1)+Validation_2(1)+Validation_3(1))/3 % 92.47%

%% sort cluster/elevation values
R1_a= cat(3, I_cluster1, RGB_map(2369:3343, 548:1758, 4)); % add elevation values to cluster image (R1)
a1= reshape(R1_a(:,:,2), 1, Size1(2)*Size1(1)); % convert elevation matrix to vector
[a1s, idx]= sort (a1, 'descend'); % sort elevation values from higher to lower and save indices
a2= reshape(R1_a(:,:,1), 1, Size1(2)*Size1(1)); % convert cluster matrix to vector
a2= a2(idx); % sort cluster values based on elevation sorted value indices
R1_as= cat(3, a2, a1s); % array with sorted cluster index-elevation values

R2_a= cat(3, I_cluster2, RGB_map(1362:2336 , 260:1470, 4)); 
a1= reshape(R2_a(:,:,2), 1, Size1(2)*Size1(1)); % convert elevation matrix to vector
[a1s, idx]= sort (a1, 'descend'); % sort elevation values from higher to lower and save indices
a2= reshape(R2_a(:,:,1), 1, Size1(2)*Size1(1)); % convert cluster matrix to vector
a2= a2(idx); % sort cluster values based on elevation sorted value indices
R2_as= cat(3, a2, a1s); % array with sorted cluster index-elevation values

R3_a= cat(3, I_cluster3, RGB_map(20:994, 372:1582, 4)); 
a1= reshape(R3_a(:,:,2), 1, Size1(2)*Size1(1)); % convert elevation matrix to vector
[a1s, idx]= sort (a1, 'descend'); % sort elevation values from higher to lower and save indices
a2= reshape(R3_a(:,:,1), 1, Size1(2)*Size1(1)); % convert cluster matrix to vector
a2= a2(idx); % sort cluster values based on elevation sorted value indices
R3_as= cat(3, a2, a1s); % array with sorted cluster index-elevation values

R4_a= cat(3, I_cluster4, RGB_map(43:1017, 2501:3711, 4)); 
a1= reshape(R4_a(:,:,2), 1, Size1(2)*Size1(1)); % convert elevation matrix to vector
[a1s, idx]= sort (a1, 'descend'); % sort elevation values from higher to lower and save indices
a2= reshape(R4_a(:,:,1), 1, Size1(2)*Size1(1)); % convert cluster matrix to vector
a2= a2(idx); % sort cluster values based on elevation sorted value indices
R4_as= cat(3, a2, a1s); % array with sorted cluster index-elevation values

R5_a= cat(3, I_cluster5, RGB_map(1362:2336, 1842:3052, 4)); 
a1= reshape(R5_a(:,:,2), 1, Size1(2)*Size1(1)); % convert elevation matrix to vector
[a1s, idx]= sort (a1, 'descend'); % sort elevation values from higher to lower and save indices
a2= reshape(R5_a(:,:,1), 1, Size1(2)*Size1(1)); % convert cluster matrix to vector
a2= a2(idx); % sort cluster values based on elevation sorted value indices
R5_as= cat(3, a2, a1s); % array with sorted cluster index-elevation values

R6_a= cat(3, I_cluster6, RGB_map(1146:2120, 3468:4678, 4)); 
a1= reshape(R6_a(:,:,2), 1, Size1(2)*Size1(1)); % convert elevation matrix to vector
[a1s, idx]= sort (a1, 'descend'); % sort elevation values from higher to lower and save indices
a2= reshape(R6_a(:,:,1), 1, Size1(2)*Size1(1)); % convert cluster matrix to vector
a2= a2(idx); % sort cluster values based on elevation sorted value indices
R6_as= cat(3, a2, a1s); % array with sorted cluster index-elevation values

%% index conversions (ATTENTION: must change each time based on k-means results)
% replace indices with random but consistent numbers (e.g. red: 5, pink: 10, black: 20)

R1_as_in= R1_as(:,:,1);
R1_as_in(R1_as_in==1)= 5;
R1_as_in(R1_as_in==2)= 10;
R1_as_in(R1_as_in==3)= 20;
R1_as= cat(3, R1_as_in, R1_as(:,:,2));
 
R2_as_in= R2_as(:,:,1);
R2_as_in(R2_as_in==1)= 10;
R2_as_in(R2_as_in==2)= 20;
R2_as_in(R2_as_in==3)= 5;
R2_as= cat(3, R2_as_in, R2_as(:,:,2));

R3_as_in= R3_as(:,:,1);
R3_as_in(R3_as_in==1)= 20;
R3_as_in(R3_as_in==2)= 5;
R3_as_in(R3_as_in==3)= 10;
R3_as= cat(3, R3_as_in, R3_as(:,:,2));
 
R4_as_in= R4_as(:,:,1);
R4_as_in(R4_as_in==1)= 10;
R4_as_in(R4_as_in==2)= 20;
R4_as_in(R4_as_in==3)= 5;
R4_as= cat(3, R4_as_in, R4_as(:,:,2));

R5_as_in= R5_as(:,:,1);
R5_as_in(R5_as_in==1)= 10;
R5_as_in(R5_as_in==2)= 20;
R5_as_in(R5_as_in==3)= 5;
R5_as= cat(3, R5_as_in, R5_as(:,:,2));

R6_as_in= R6_as(:,:,1);
R6_as_in(R6_as_in==1)= 5;
R6_as_in(R6_as_in==2)= 20;
R6_as_in(R6_as_in==3)= 10;
R6_as= cat(3, R6_as_in, R6_as(:,:,2));

 
%% elevation ranges (ATTENTION: change indices each time)
% cluster 1
R1_1= find(R1_as(:,:,1)==20); % pixel indices belonging to cluster 1
Elev= R1_as(:,:,2); % isolate elevation vector
R1_1_elev=  Elev(R1_1); % elevation values for cluster 1
csvwrite('output_path', R1_1_elev)
max(R1_1_elev);
min(R1_1_elev);
mean(R1_1_elev);
std(R1_1_elev);
hist(R1_1_elev, 30); % 30 bins

R1_2= find(R1_as(:,:,1)==10); % cluster 2
R1_2_elev=  Elev(R1_2); 
csvwrite('output_path', R1_2_elev)
max(R1_2_elev);
min(R1_2_elev);
mean(R1_2_elev);
std(R1_2_elev);
hist(R1_2_elev, 30);

R1_3= find(R1_as(:,:,1)==5); % cluster 3
R1_3_elev=  Elev(R1_3); 
csvwrite('output_path', R1_3_elev)
max(R1_3_elev);
min(R1_3_elev);
mean(R1_3_elev);
std(R1_3_elev);
hist(R1_3_elev, 30);

% cluster 2
R2_1= find(R2_as(:,:,1)==10); % pixel indices belonging to cluster 1
Elev= R2_as(:,:,2); % isolate elevation vector
R2_1_elev=  Elev(R2_1); % elevation values for cluster 1
csvwrite('output_path', R2_1_elev)
max(R2_1_elev);
min(R2_1_elev);
mean(R2_1_elev);
std(R2_1_elev);
hist(R2_1_elev, 30);

R2_2= find(R2_as(:,:,1)==20); % cluster 2
R2_2_elev=  Elev(R2_2); 
csvwrite('output_path', R2_2_elev)
max(R2_2_elev);
min(R2_2_elev);
mean(R2_2_elev);
std(R2_2_elev);
hist(R2_2_elev, 30);

R2_3= find(R2_as(:,:,1)==5); % cluster 3
R2_3_elev=  Elev(R2_3); 
csvwrite('output_path', R2_3_elev)
max(R2_3_elev);
min(R2_3_elev);
mean(R2_3_elev);
std(R2_3_elev);
hist(R2_3_elev, 30);

% cluster 3
R3_1= find(R3_as(:,:,1)==5); % pixel indices belonging to cluster 1
Elev= R3_as(:,:,2); % isolate elevation vector
R3_1_elev=  Elev(R3_1); % elevation values for cluster 1
csvwrite('output_path', R3_1_elev)
max(R3_1_elev);
min(R3_1_elev);
mean(R3_1_elev);
std(R3_1_elev);
hist(R3_1_elev, 30);

R3_2= find(R3_as(:,:,1)==10); % cluster 2
R3_2_elev=  Elev(R3_2); 
csvwrite('output_path', R3_2_elev)
max(R3_2_elev);
min(R3_2_elev);
mean(R3_2_elev);
std(R3_2_elev);
hist(R3_2_elev, 30);

R3_3= find(R3_as(:,:,1)==20); % cluster 3
R3_3_elev=  Elev(R3_3); 
csvwrite('output_path', R3_3_elev)
max(R3_3_elev);
min(R3_3_elev);
mean(R3_3_elev);
std(R3_3_elev);
hist(R3_3_elev, 30);

% cluster 4
R4_1= find(R4_as(:,:,1)==10); % pixel indices belonging to cluster 1
Elev= R4_as(:,:,2); % isolate elevation vector
R4_1_elev=  Elev(R4_1); % elevation values for cluster 1
csvwrite('output_path', R4_1_elev)
max(R4_1_elev);
min(R4_1_elev);
mean(R4_1_elev);
std(R4_1_elev);
hist(R4_1_elev, 30);

R4_2= find(R4_as(:,:,1)==5); % cluster 2
R4_2_elev=  Elev(R4_2);
csvwrite('output_path', R4_2_elev)
max(R4_2_elev);
min(R4_2_elev);
mean(R4_2_elev);
std(R4_2_elev);
hist(R4_2_elev, 30);

R4_3= find(R4_as(:,:,1)==20); % cluster 3
R4_3_elev=  Elev(R4_3); 
csvwrite('output_path', R4_3_elev)
max(R4_3_elev);
min(R4_3_elev);
mean(R4_3_elev);
std(R4_3_elev);
hist(R4_3_elev, 30);

% cluster 5
R5_1= find(R5_as(:,:,1)==10); % pixel indices belonging to cluster 1
Elev= R5_as(:,:,2); % isolate elevation vector
R5_1_elev=  Elev(R5_1); % elevation values for cluster 1
csvwrite('output_path', R5_1_elev)
max(R5_1_elev);
min(R5_1_elev);
mean(R5_1_elev);
std(R5_1_elev);
hist(R5_1_elev, 30);

R5_2= find(R5_as(:,:,1)==20); % cluster 2
R5_2_elev=  Elev(R5_2); 
csvwrite('output_path', R5_2_elev)
max(R5_2_elev);
min(R5_2_elev);
mean(R5_2_elev);
std(R5_2_elev);
hist(R5_2_elev, 30);

R5_3= find(R5_as(:,:,1)==5); % cluster 3
R5_3_elev=  Elev(R5_3); 
csvwrite('output_path', R5_3_elev)
max(R5_3_elev);
min(R5_3_elev);
mean(R5_3_elev);
std(R5_3_elev);
hist(R5_3_elev, 30);

% cluster 6
R6_1= find(R6_as(:,:,1)==20); % pixel indices belonging to cluster 1
Elev= R6_as(:,:,2); % isolate elevation vector
R6_1_elev=  Elev(R6_1); % elevation values for cluster 1
csvwrite('output_path', R6_1_elev)
max(R6_1_elev);
min(R6_1_elev);
mean(R6_1_elev);
std(R6_1_elev);
hist(R6_1_elev, 30);

R6_2= find(R6_as(:,:,1)==5); % cluster 2
R6_2_elev=  Elev(R6_2);
csvwrite('output_path', R6_2_elev)
max(R6_2_elev);
min(R6_2_elev);
mean(R6_2_elev);
std(R6_2_elev);
hist(R6_2_elev, 30);

R6_3= find(R6_as(:,:,1)==10); % cluster 3
R6_3_elev=  Elev(R6_3); 
csvwrite('output_path', R6_3_elev)
max(R6_3_elev);
min(R6_3_elev);
mean(R6_3_elev);
std(R6_3_elev);
hist(R6_3_elev, 30);

% all images (ATTENTION: each image may have different indices per species)
All= cat(1, R1_as, R2_as, R3_as, R4_as, R5_as, R6_as); % unite all images by rows
All= reshape(All, 1, 6*1180725, 2); % reshape to array of two vectors
ElevAll= All(:,:,2);
ClusterAll= All(:,:,1);

All_1= find(All(:,:,1)==5); % cluster 1
All_1_elev=  ElevAll(All_1);
csvwrite('output_path', All_1_elev)
max(All_1_elev);
min(All_1_elev);
mean(All_1_elev);
std(All_1_elev);
hist(All_1_elev, 30); % average elevation distribution of cluster 1 per all images

All_2= find(All(:,:,1)==10); % cluster 2
All_2_elev=  ElevAll(All_2);
csvwrite('output_path', All_2_elev)
max(All_2_elev);
min(All_2_elev);
mean(All_2_elev);
std(All_2_elev);
hist(All_2_elev, 30);

All_3= find(All(:,:,1)==20); % cluster 3
All_3_elev=  ElevAll(All_3);
csvwrite('output_path', All_3_elev)
max(All_3_elev);
min(All_3_elev);
mean(All_3_elev);
std(All_3_elev);
hist(All_3_elev, 30);
%% elevation range proportions

Gradient= transpose(0.05:0.05:1); % turn rows into columns to concatenate later
Elev1= RGB_map(2369:3343, 548:1758, 4);
Elev1_prop= [];
Elev1_abund=[];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(Elev1>=i & Elev1<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(a)/numel(Elev1); % proportion of pixels belonging to i and i+0.05
    b1= numel(a)
    Elev1_prop= cat(1, Elev1_prop, a1); % add all proportion values (from lowest to highest) in vector
    Elev1_abund= cat(1, Elev1_abund, b1);
end
Elev1_prop= cat(2, Elev1_prop, Gradient);
Elev1_abund= cat(2, Elev1_abund, Gradient);
csvwrite('output_path', Elev1_prop)
csvwrite('output_path', Elev1_abund)

Elev2= RGB_map(1362:2336, 548:1758, 4);
Elev2_prop= [];
Elev2_abund=[];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(Elev2>=i & Elev2<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(a)/numel(Elev2); % proportion of pixels belonging to i and i+0.05
    b1= numel(a)
    Elev2_prop= cat(1, Elev2_prop, a1); % add all proportion values (from lowest to highest) in vector
    Elev2_abund= cat(1, Elev2_abund, b1);

end
Elev2_prop= cat(2, Elev2_prop, Gradient);
Elev2_abund= cat(2, Elev1_abund, Gradient);
csvwrite('output_path', Elev2_prop)
csvwrite('output_path', Elev2_abund)


Elev3= RGB_map(20:994, 548:1758, 4);
Elev3_prop= [];
Elev3_abund=[];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(Elev3>=i & Elev3<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(a)/numel(Elev3); % proportion of pixels belonging to i and i+0.05
    b1= numel(a)
    Elev3_prop= cat(1, Elev3_prop, a1); % add all proportion values (from lowest to highest) in vector
    Elev3_abund= cat(1, Elev3_abund, b1);
end
Elev3_prop= cat(2, Elev3_prop, Gradient);
Elev3_abund= cat(2, Elev3_abund, Gradient);
csvwrite('output_path', Elev3_prop)
csvwrite('output_path', Elev3_abund)


Elev4= RGB_map(43:1017, 548:1758, 4);
Elev4_prop= [];
Elev4_abund=[];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(Elev4>=i & Elev4<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(a)/numel(Elev4); % proportion of pixels belonging to i and i+0.05
        b1= numel(a)
    Elev4_prop= cat(1, Elev4_prop, a1); % add all proportion values (from lowest to highest) in vector
    Elev4_abund= cat(1, Elev4_abund, b1);
end
Elev4_prop= cat(2, Elev4_prop, Gradient);
Elev4_abund= cat(2, Elev4_abund, Gradient);
csvwrite('output_path', Elev4_prop)
csvwrite('output_path', Elev4_abund)


Elev5= RGB_map(1362:2336, 548:1758, 4);
Elev5_prop= [];
Elev5_abund=[];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(Elev5>=i & Elev5<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(a)/numel(Elev5); % proportion of pixels belonging to i and i+0.05
    b1= numel(a)
    Elev5_prop= cat(1, Elev5_prop, a1); % add all proportion values (from lowest to highest) in vector
    Elev5_abund= cat(1, Elev5_abund, b1);
end
Elev5_prop= cat(2, Elev5_prop, Gradient);
Elev5_abund= cat(2, Elev5_abund, Gradient);
csvwrite('output_path', Elev5_abund)
csvwrite('output_path', Elev5_prop)


Elev6= RGB_map(1146:2120, 548:1758, 4);
Elev6_prop= [];
Elev6_abund= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(Elev6>=i & Elev6<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(a)/numel(Elev6); % proportion of pixels belonging to i and i+0.05
    b1= numel(a)
    Elev6_prop= cat(1, Elev6_prop, a1); % add all proportion values (from lowest to highest) in vector
    Elev6_abund= cat(1, Elev6_abund, b1);
end
Elev6_prop= cat(2, Elev6_prop, Gradient);
Elev6_abund= cat(2, Elev6_abund, Gradient);
csvwrite('output_path', Elev6_abund)
csvwrite('output_path', Elev6_prop)

%% rate of relative abundance change per unit of relative altitude
% image 1
Gradient= transpose(0.05:0.05:1); % turn rows into columns to concatenate later
R1_gradient_1= [];
R1_gradient_2= [];
R1_gradient_3= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(R1_as(:,:,2)>=i & R1_as(:,:,2)<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(find(R1_as(a)==20))/numel(a); % proportion of pixels belonging to cluster 1 between i and i+1
    a2= numel(find(R1_as(a)==10))/numel(a);
    a3= numel(find(R1_as(a)==5))/numel(a);
    R1_gradient_1= cat(1, R1_gradient_1, a1); % add all proportion values (from lowest to highest) in vector
    R1_gradient_2= cat(1, R1_gradient_2, a2);
    R1_gradient_3= cat(1, R1_gradient_3, a3);
end
R1_gradient_1= cat(2, R1_gradient_1, Gradient); % left column corresponds to relative abundance and right to elevation range
csvwrite('output_path', R1_gradient_1)
plot(R1_gradient_1(:,2), R1_gradient_1(:,1))
R1_gradient_2= cat(2, R1_gradient_2, Gradient);
csvwrite('output_path', R1_gradient_2)
plot(R1_gradient_2(:,2), R1_gradient_2(:,1))
R1_gradient_3= cat(2, R1_gradient_3, Gradient);
csvwrite('output_path', R1_gradient_3)
plot(R1_gradient_3(:,2), R1_gradient_3(:,1))

% image 2
R2_gradient_1= [];
R2_gradient_2= [];
R2_gradient_3= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(R2_as(:,:,2)>=i & R2_as(:,:,2)<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(find(R2_as(a)==10))/numel(a); % proportion of pixels belonging to cluster 1 between i and i+1
    a2= numel(find(R2_as(a)==20))/numel(a);
    a3= numel(find(R2_as(a)==5))/numel(a);
    R2_gradient_1= cat(1, R2_gradient_1, a1); % add all proportion values (from lowest to highest) in vector
    R2_gradient_2= cat(1, R2_gradient_2, a2);
    R2_gradient_3= cat(1, R2_gradient_3, a3);
end
R2_gradient_1= cat(2, R2_gradient_1, Gradient); % left column corresponds to relative abundance and right to elevation range
csvwrite('output_path', R2_gradient_1)
plot(R2_gradient_1(:,2), R2_gradient_1(:,1))
R2_gradient_2= cat(2, R2_gradient_2, Gradient);
csvwrite('output_path', R2_gradient_2)
plot(R2_gradient_2(:,2), R2_gradient_2(:,1))
R2_gradient_3= cat(2, R2_gradient_3, Gradient);
csvwrite('output_path', R2_gradient_3)
plot(R2_gradient_3(:,2), R2_gradient_3(:,1))

% image 3

R3_gradient_1= [];
R3_gradient_2= [];
R3_gradient_3= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(R3_as(:,:,2)>=i & R3_as(:,:,2)<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(find(R3_as(a)==5))/numel(a); % proportion of pixels belonging to cluster 1 between i and i+1
    a2= numel(find(R3_as(a)==10))/numel(a);
    a3= numel(find(R3_as(a)==20))/numel(a);
    R3_gradient_1= cat(1, R3_gradient_1, a1); % add all proportion values (from lowest to highest) in vector
    R3_gradient_2= cat(1, R3_gradient_2, a2);
    R3_gradient_3= cat(1, R3_gradient_3, a3);
end
R3_gradient_1= cat(2, R3_gradient_1, Gradient); % left column corresponds to relative abundance and right to elevation range
csvwrite('output_path', R3_gradient_1)
plot(R2_gradient_1(:,2), R3_gradient_1(:,1))
R3_gradient_2= cat(2, R3_gradient_2, Gradient);
csvwrite('output_path', R3_gradient_2)
plot(R2_gradient_2(:,2), R3_gradient_2(:,1))
R3_gradient_3= cat(2, R3_gradient_3, Gradient);
csvwrite('output_path', R3_gradient_3)
plot(R3_gradient_3(:,2), R3_gradient_3(:,1))

% image 4

R4_gradient_1= [];
R4_gradient_2= [];
R4_gradient_3= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(R4_as(:,:,2)>=i & R4_as(:,:,2)<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(find(R4_as(a)==10))/numel(a); % proportion of pixels belonging to cluster 1 between i and i+1
    a2= numel(find(R4_as(a)==5))/numel(a);
    a3= numel(find(R4_as(a)==20))/numel(a);
    R4_gradient_1= cat(1, R4_gradient_1, a1); % add all proportion values (from lowest to highest) in vector
    R4_gradient_2= cat(1, R4_gradient_2, a2);
    R4_gradient_3= cat(1, R4_gradient_3, a3);
end
R4_gradient_1= cat(2, R4_gradient_1, Gradient); % left column corresponds to relative abundance and right to elevation range
csvwrite('output_path', R4_gradient_1)
plot(R2_gradient_1(:,2), R4_gradient_1(:,1))
R4_gradient_2= cat(2, R4_gradient_2, Gradient);
csvwrite('output_path', R4_gradient_2)
plot(R4_gradient_2(:,2), R4_gradient_2(:,1))
R4_gradient_3= cat(2, R4_gradient_3, Gradient);
csvwrite('output_path', R4_gradient_3)
plot(R4_gradient_3(:,2), R4_gradient_3(:,1))

% image 5

R5_gradient_1= [];
R5_gradient_2= [];
R5_gradient_3= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(R5_as(:,:,2)>=i & R5_as(:,:,2)<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(find(R5_as(a)==10))/numel(a); % proportion of pixels belonging to cluster 1 between i and i+1
    a2= numel(find(R5_as(a)==20))/numel(a);
    a3= numel(find(R5_as(a)==5))/numel(a);
    R5_gradient_1= cat(1, R5_gradient_1, a1); % add all proportion values (from lowest to highest) in vector
    R5_gradient_2= cat(1, R5_gradient_2, a2);
    R5_gradient_3= cat(1, R5_gradient_3, a3);
end
R5_gradient_1= cat(2, R5_gradient_1, Gradient); % left column corresponds to relative abundance and right to elevation range
csvwrite('output_path', R5_gradient_1)
plot(R2_gradient_1(:,2), R5_gradient_1(:,1))
R5_gradient_2= cat(2, R5_gradient_2, Gradient);
csvwrite('output_path', R5_gradient_2)
plot(R5_gradient_2(:,2), R5_gradient_2(:,1))
R5_gradient_3= cat(2, R5_gradient_3, Gradient);
csvwrite('output_path', R5_gradient_3)
plot(R5_gradient_3(:,2), R5_gradient_3(:,1))

% image 6

R6_gradient_1= [];
R6_gradient_2= [];
R6_gradient_3= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(R6_as(:,:,2)>=i & R6_as(:,:,2)<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(find(R6_as(a)==20))/numel(a); % proportion of pixels belonging to cluster 1 between i and i+1
    a2= numel(find(R6_as(a)==5))/numel(a);
    a3= numel(find(R6_as(a)==10))/numel(a);
    R6_gradient_1= cat(1, R6_gradient_1, a1); % add all proportion values (from lowest to highest) in vector
    R6_gradient_2= cat(1, R6_gradient_2, a2);
    R6_gradient_3= cat(1, R6_gradient_3, a3);
end
R6_gradient_1= cat(2, R6_gradient_1, Gradient); % left column corresponds to relative abundance and right to elevation range
csvwrite('output_path', R6_gradient_1)
plot(R6_gradient_1(:,2), R6_gradient_1(:,1))
R6_gradient_2= cat(2, R6_gradient_2, Gradient);
csvwrite('output_path', R6_gradient_2)
plot(R6_gradient_2(:,2), R6_gradient_2(:,1))
R6_gradient_3= cat(2, R6_gradient_3, Gradient);
csvwrite('output_path', R6_gradient_3)
plot(R6_gradient_3(:,2), R6_gradient_3(:,1))

% all clusters (ATTENTION: each image may have different indices per species)

All_gradient_1= [];
All_gradient_2= [];
All_gradient_3= [];
for i=0:0.05:0.95 % break elevation range to 5% values
    a= find(ElevAll>=i & ElevAll<=i+0.05); % find elevation values from i to i+0.05 (0.05 difference)
    a1= numel(find(ClusterAll(a)==5))/numel(a); % proportion of pixels belonging to cluster 1 between i and i+1
    a2= numel(find(ClusterAll(a)==10))/numel(a);
    a3= numel(find(ClusterAll(a)==20))/numel(a);
    All_gradient_1= cat(1, All_gradient_1, a1); % add all proportion values (from lowest to highest) in vector
    All_gradient_2= cat(1, All_gradient_2, a2);
    All_gradient_3= cat(1, All_gradient_3, a3);
end
All_gradient_1= cat(2, All_gradient_1, Gradient); % left column corresponds to relative abundance and right to elevation range
csvwrite('output_path', All_gradient_1)
plot(All_gradient_1(:,2), All_gradient_1(:,1))
All_gradient_2= cat(2, All_gradient_2, Gradient);
csvwrite('output_path', All_gradient_2)
plot(All_gradient_2(:,2), All_gradient_2(:,1))
All_gradient_3= cat(2, All_gradient_3, Gradient);
csvwrite('output_path', All_gradient_3)
plot(All_gradient_3(:,2), All_gradient_3(:,1))

All_gradient_123= cat(2, All_gradient_1(:,1),All_gradient_2(:,1),All_gradient_3(:,1), Gradient);
csvwrite('output_path', All_gradient_123)
