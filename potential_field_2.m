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
w = 1/50;  %Wider contours avoid sharp changes in projection direction for points at the intersections, 1/50 works well

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
axis([-50 200 -50 150 -10 10])
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



%{
coeffs_x = polyfit(X(1,1:10),1:10,1);
coeffs_y = polyfit(Y(1:10,1)',1:10,1);
ind_x = int16((coeffs_x(1)*(raw_pos(i,1)))+coeffs_x(2));
ind_y = int16((coeffs_y(1)*(raw_pos(i,2)))+coeffs_y(2));
function index = x_ind(x,coeffs_x)
index = (coeffs_x(1)*(x))+coeffs_x(2);
end
ind_x
%}

tic;
f_grad = [];
for i=1:length(raw_pos(:,1))
    t = (X(:)==raw_pos(i,1))&(Y(:)==raw_pos(i,2));
    indt = find(t);
    f_grad(i,:) = [px(indt) py(indt)];
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
new_pos = [];
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
    
     %{
    t_points = segmentCoords([1 2 4],3:4);
    if    ~((all_dir_pos(index,1)>segs_sorted(index,1))...
                   &(all_dir_pos(index,1)<segs_sorted(index,2))...
                   &(all_dir_pos(index,2)>segs_sorted(index,3))...
                   &(all_dir_pos(index,2)<segs_sorted(index,4)))
               
               d_from_t = all_dir_pos(index,1:2)-t_points;
               distance = [];
               for k = 1:size(d_from_t,1)
                   distance(k) = norm(d_from_t(k,:));
               end
               [val closest_t] = min(distance);
               new_pos(i,:) = t_points(closest_t,:);
               
    else             
                new_pos(i,:) = all_dir_pos(index,1:2);
    end
     %}
    
        new_pos(i,:) = all_dir_pos(index,1:2);
    if ~(mod(i,1000))
        %display percent progress
        disp([num2str(round((i/poslength)*100)),'%']);
    end
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