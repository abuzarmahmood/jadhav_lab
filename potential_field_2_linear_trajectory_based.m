%% Potential and Gradient
clear all;
load('wenbos_variables.mat');
%Clean start

seg_start = segmentCoords(:,1:2);
seg_end = segmentCoords(:,3:4);

segmentCoords = round(segmentCoords,1);
seg_dirs = segmentCoords(:,3:4) - segmentCoords(:,1:2);
for i=1:size(seg_dirs,1)
    seg_dirs(i,:) = seg_dirs(i,:)/norm(seg_dirs(i,:));
end

load('KL8trajinfo01.mat');
traj_wells =  trajinfo{1, 1}{1, 10}.wellstend  ;
traj_time = trajinfo{1, 1}{1, 10}.trajtime;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%Generating the potential
%%%%%%%%%%%%%%%%%%%%%%%%%

raw_pos= round([pos(:,2) pos(:,3)],1);
buffer = 20; % how much extra around edges
den = 0.1;

[X,Y] = meshgrid((min(raw_pos(:,1))-buffer):den:(max(raw_pos(:,1))+buffer),(min(raw_pos(:,2))-buffer):den:(max(raw_pos(:,2))+buffer));  % covers all points with grid spacing = den
X = round(X,1);
Y = round(Y,1);


w = 1/100;

f_pots = zeros(size(X,1),size(X,2),2); % for 2 trajectories
tic;
for path = 1:2
        trajs = [1 2 3; 1 4 5];
        coords_used = segmentCoords(trajs(path,:),:);
        pots = zeros(size(X,1),size(X,2),3); % for 3 segments
    
               
    for i = 1:size(coords_used,1)
        if ~mod(i,2) % even segment = horizontal
        coeffs = polyfit(coords_used(i,[1,3]),coords_used(i,[2,4]),1); % fit for y = mx+c
        pots(:,:,i) = (((1+exp(-1.*((((coeffs(1).*X)+coeffs(2))-Y).^2).*w)).^-1)-1);
        %pots(
        else         % odd segment = vertical
        coeffs = polyfit(coords_used(i,[2,4]),coords_used(i,[1,3]),1);  % fit for x = my+c
        pots(:,:,i) = (((1+exp(-1.*((((coeffs(1).*Y)+coeffs(2))-X).^2).*w)).^-1)-1);        
        end
    end

    %###Wells at intersections need to be addressed###%
    sum_pots = sum(pots,3);
    %sum_pots(sum_pots < median(all_min)) = median(all_min);
    f_pots(:,:,path) = sum_pots;
    
    [px,py] = gradient(f_pots(:,:,path),0.1); % Calculate gradient
    grad{path,1} = px;
    grad{path,2} = py;
end

%{
figure
contour(X,Y,f_pots)
s = surf(X,Y,f_pots(:,:,1));
%axis([-50 200 -50 150 -10 10])
set(s,'LineStyle','none')
hold on

quiver(X,Y,px,py);
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

%% 
%%%%%%%%%%%%%
%%Projection
%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Gradient lookup for every point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Messing with seg_dirs and seg_start
if length(segmentLength)==5
    seg_dirs(2,:) = (seg_dirs(2,:)-seg_dirs(4,:))./2;
    %seg_dirs(4,1) = seg_dirs(4,1)*(-1);
    seg_dirs = seg_dirs([1 2 3 5],:);
    seg_start = seg_start([1 2 3 5],:);
    segmentLength = segmentLength([1 2 3 5]);
end
%}

%Ideally, divide raw_pos into trajectory-based groups and then find f_grad
%for each group separately.
tic; 
f_grad = {};
for i=1:length(raw_pos(:,1))
    for j = 1:size(grad,1)
        t = (X(:)==raw_pos(i,1))&(Y(:)==raw_pos(i,2));
        indt = find(t);
        f_grad{j}(i,:) = [grad{j,1}(indt) grad{j,2}(indt)];
        f_grad{j}(i,:) = f_grad{j}(i,:)/norm(f_grad{j}(i,:));

        if ~(mod(i,1000))
            %display percent progress
            disp([num2str(round((i/poslength)*100)),'%']);
        end
    end
end
toc;

segs_sorted = [sort(segmentCoords(:,[1 3]),2) sort(segmentCoords(:,[2 4]),2)];
segs_sorted([1 3 5],4) = max(raw_pos(:,1)+2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Dividing position and gradient by trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

trajtimes_1to2 = traj_time(find(sum((traj_wells == [1 2])| (traj_wells == [2 1]),2)-1),:);
trajtimes_1to3 = traj_time(find(sum((traj_wells == [1 3])| (traj_wells == [3 1]),2)-1),:);
trajtimes = {trajtimes_1to2 trajtimes_1to3};

traj_pos = {[],[]};
quiv = {[],[]};

%Calculate raw position and gradient vector for points in both trajectory
%groups.
for traj = 1:2
    for i = 1:size(trajtimes{traj},1)
        start_time = trajtimes{traj}(i,1);
        end_time = trajtimes{traj}(i,2);
        start_ind = find(pos(:,1)==start_time);
        end_ind = find(pos(:,1)==end_time);
        
        traj_pos{traj} = [traj_pos{traj}; raw_pos(start_ind:end_ind,:),linspace(start_time,end_time,end_ind-start_ind+1)'];
        quiv{traj} = [quiv{traj};f_grad{traj}(start_ind:end_ind,:),linspace(start_time,end_time,end_ind-start_ind+1)'];
        
    end
end

%%
%%%%%%%%%%%%%%%%%%
%%Final Projection
%%%%%%%%%%%%%%%%%%

new_pos = {[],[]};
tic;
    for lr = 1:size(trajs,1)   %Loop over both U trajectories (left or right)

        for i=1:size(traj_pos{lr},1) % Loop through all points in each U trajectory
            
            all_dir_pos = [];
            proj = [];
            
            for j = 1:size(seg_dirs,1) % Loop through all segments
                
                M = [seg_dirs(j,:)',quiv{lr}(i,1:2)'];
                trans_mat = pinv(M'*M)*M';
                p_rel = traj_pos{lr}(i,1:2)'-seg_start(j,:)';
                proj(j,:) = trans_mat*p_rel;

                all_dir_pos(j,1:2) = seg_start(j,:) + proj(j,1).*seg_dirs(j,:);
                all_dir_pos(j,3) = norm(traj_pos{lr}(i,1:2) - all_dir_pos(j,1:2));
            end

            %To make sure the projection is inside the particular segment
             [val index] = min(all_dir_pos(:,3));


                new_pos{lr} = [new_pos{lr};all_dir_pos(index,1:2)];
            if ~(mod(i,1000))
                %display percent progress
                disp([num2str(round((i/poslength)*100)),'%']);
            end
            
        end
    end
toc;
    
%{
level = 2; %1 or 2
scatter(traj_pos{level}(:,1),traj_pos{level}(:,2),2)
hold on
plot(new_pos{level}(:,1),new_pos{level}(:,2))
contour(X,Y,f_pots(:,:,level))
q = quiver(traj_pos{level}(:,1),traj_pos{level}(:,2),-quiv{level}(:,1),-quiv{level}(:,2));
q.AutoScaleFactor = 0.3;
%}