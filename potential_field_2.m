%% Potential and Gradient
%clear all;
%load('wenbos_variables.mat');
%Clean start
seg_start = segmentCoords(:,1:2);
seg_end = segmentCoords(:,3:4);
seg_length = segmentLength;

segmentCoords = round(segmentCoords,1);
seg_dirs = segmentCoords(:,3:4) - segmentCoords(:,1:2);
for i=1:size(seg_dirs,1)
    seg_dirs(i,:) = seg_dirs(i,:)/norm(seg_dirs(i,:));
    
end
%%%%%%%%%%%%%%%%%
    seg_start2 = segmentCoords(:,1:2);
    seg_end2 = segmentCoords(:,3:4);
    seg_length2 = segmentLength;

    segmentCoords2 = round(segmentCoords,1);
    seg_dirs2 = segmentCoords(:,3:4) - segmentCoords(:,1:2);
    for i=1:size(seg_dirs2,1)
        seg_dirs2(i,:) = seg_dirs2(i,:)/norm(seg_dirs2(i,:));
    end
%%%%%%%%%%%%%%%%%

%Messing with seg_dirs and seg_start
if size(seg_dirs,1)==5
    seg_dirs(2,:) = (seg_dirs(2,:)-seg_dirs(4,:))./2;
    %seg_dirs(4,1) = seg_dirs(4,1)*(-1);
    seg_dirs = seg_dirs([1 2 3 5],:);
    seg_start = seg_start([1 2 3 5],:);
    seg_length = seg_length([1 2 3 5]);
end
%}


%%

raw_pos= round([pos(:,2) pos(:,3)],1);
buffer = 20; % how much extra around edges
den = 0.1;

[X,Y] = meshgrid((min(raw_pos(:,1))-buffer):den:(max(raw_pos(:,1))+buffer),(min(raw_pos(:,2))-buffer):den:(max(raw_pos(:,2))+buffer));  % covers all points with grid spacing = den
X = round(X,1);
Y = round(Y,1);
pots = zeros(size(X,1),size(X,2),5); % for 5 segments
w = 1/100;  %Wider contours avoid sharp changes in projection direction for points at the intersections, 1/50 works well

for i = 1:length(coord)
    if ~mod(i,2) % even segment = horizontal
        coeffs = polyfit(coord{i}(:,1),coord{i}(:,2),1); % fit for y = mx+c
        pots(:,:,i) = (((1+exp(-1.*((((coeffs(1).*X)+coeffs(2))-Y).^2).*w)).^-1)-1);
        %pots(
    else         % odd segment = vertical
        coeffs = polyfit(coord{i}(:,2),coord{i}(:,1),1);  % fit for x = my+c
        pots(:,:,i) = (((1+exp(-1.*((((coeffs(1).*Y)+coeffs(2))-X).^2).*w)).^-1)-1);        
    end
end

f_pots = sum(pots,3);
[px,py] = gradient(f_pots); % Calculate gradient

%{
figure
contour(X,Y,f_pots)
s = surf(X,Y,f_pots)
%axis([-50 200 -50 150 -10 10])
set(s,'LineStyle','none')
hold on

quiver(X,Y,log(abs(px)),log(abs(py)));
scatter(raw_pos(:,1),raw_pos(:,2),5)
hold off
%}

%{
figure
contour(X,Y,f_pots)
hold on
scatter(raw_pos(i,1),raw_pos(i,2));
q = quiver(raw_pos(i,1),raw_pos(i,2),-f_grad(i,1),-f_grad(i,2));
scatter(new_pos(i,1),new_pos(i,2));
hold off
%s = q.AutoScaleFactor;
%q.AutoScaleFactor = 1000;
%}

%% Projection

%% 


tic;

x_fit = polyfit(X(1,1:10),1:10,1);
y_fit = polyfit(Y(1:10,1)',1:10,1);

f_grad = [];
for i=1:length(raw_pos(:,1))
    x_val = int16(x_fit(1)*raw_pos(i,1) + x_fit(2)); %x index gives columns and y index gives rows
    y_val = int16(y_fit(1)*raw_pos(i,2) + y_fit(2));
    
    f_grad(i,:) = [px(y_val,x_val) py(y_val,x_val)];
    f_grad(i,:) = f_grad(i,:)/norm(f_grad(i,:));
    
    if ~(mod(i,1000))
        %display percent progress
        disp([num2str(round((i/poslength)*100)),'%']);
    end
end
toc;

segs_sorted = [sort(segmentCoords(:,[1 3]),2) sort(segmentCoords(:,[2 4]),2)];
segs_sorted([1 3 5],4) = max(raw_pos(:,1)+2);

%%
tic;
new_pos = zeros(size(raw_pos,1),3);
for i=1:length(raw_pos(:,1))
    all_dir_pos = [];
    proj = [];
    for j = 1:size(seg_dirs,1)
        M = [seg_dirs(j,:)',f_grad(i,:)'];
        trans_mat = pinv(M'*M)*M';
        p_rel = raw_pos(i,:)'-seg_start(j,:)';
        proj(j,:) = trans_mat*p_rel;
         
        all_dir_pos(j,1:2) = seg_start(j,:) + proj(j,1).*seg_dirs(j,:);
        all_dir_pos(j,3) = norm(raw_pos(i,:) - all_dir_pos(j,1:2));
    end
    
    %To make sure the projection is inside the particular segment
     [val index] = min(all_dir_pos(:,3)); 
     
     
    
        new_pos(i,1:2) = all_dir_pos(index,1:2);
        
        if index==2 %%To allow use of old distanceTable
         if new_pos(i,1)<=seg_end(1,1)
             index = 2;
         else 
             index = 5;
         end
        end
     
        new_pos(i,3) = index;  %%for calculation of lindist
        
        
    if ~(mod(i,1000))
        %display percent progress
        disp([num2str(round((i/poslength)*100)),'%']);
    end
end
toc;

%%
%%%%%%%%%%%%%%%%%%%%
%lindist calculation
%%%%%%%%%%%%%%%%%%%%
%1)Index of segment projection is on
%2)Is segment start point forwards or backwards of all of the 3 wells
%3)Add or subtact distnace of projection from segment start depending on
%forward or backwards

% **Have to calculate disance from resepective starting point since I can't
% use proj because of splitting horizontal segment.

%New segment layout
%   *    *    *
%   |    |    |
%   3    1    4
%   |    |    |
%   |    |    |
%   *-2--*-5--*
%
tic;

new_dist_table = distanceTable(:,[1 2 3 5 4]);
new_seg_start = seg_start2([1 2 3 5 4],:);
new_seg_length = seg_length2([1 2 3 5 4]);
lindist_abu = zeros(size(new_pos,1),3);
%segmentDirection is same since flipping 4 and 5 has no effect. 

for i = 1:size(new_pos,1)
    index = new_pos(i,3);
    dist_from_start = norm(new_pos(i,1:2)-new_seg_start(index,:));
    distToSegment = new_dist_table(:,index)';
    for wellcount = 1:size(new_dist_table,1)
        
        
        
                    if  (segmentDirection(wellcount,index) == 1)
                        %it is aligned in the foreward direction
                        segmentdist(1,wellcount)  = dist_from_start;
                    else
                        %it is aligned in the backward direction
                        segmentdist(1,wellcount)  = seg_length2(index)-dist_from_start;
                    end
    end
    lindist_abu(i,1:wellcount) = distToSegment + segmentdist;
end

toc;
%{
figure
contour(X,Y,f_pots)
hold on
scatter(raw_pos(:,1),raw_pos(:,2));
q = quiver(raw_pos(:,1),raw_pos(:,2),-f_grad(:,1),-f_grad(:,2));
plot(new_pos(:,1),new_pos(:,2));
hold off
%s = q.AutoScaleFactor;
%q.AutoScaleFactor = 1000;
%}