%
%%
%Setup: Code imported from Wenbo/Mathias (till potential and gradient)
%Need to load files pos amd task, setup segmentCoords and distanceTable.
    tic;
    clear all;
    close all;
    clc;
    %%
    % add the below folders to your MATLAB path
    addpath(genpath('C:\Users\Justin-Admin\Downloads\Linpos code\Src_Matlab\'));
    addpath(genpath('C:\Users\Justin-Admin\Downloads\Linpos code\usrlocal\'));

    %%
    %choose your animal
    animalprefix = 'KL8';
    day = 1;
    ep = 10;

    % change the directory of the animal folder accordingly
    if strcmp(animalprefix,'ER1')
        dir = 'D:\SingledayExp\ER1_NEW_direct2\';
    elseif strcmp(animalprefix,'KL8')
        dir = 'C:\Users\Justin-Admin\Downloads\Linpos code\';
    end

    pos = loaddatastruct(dir, animalprefix, 'pos', day); % get pos file
    task = loaddatastruct(dir, animalprefix, 'task', day); % get task info

    %%
    % parameters for linearization
    maxsegdiff = 100;
    maxv = 300;
    %%
    % make sure we have the direction information
    toknum = isdatafield(pos{day}{ep}.fields, 'dir');
    if (~toknum)
        Error('No direction field in pos');
    end

    pos = pos{day}{ep}.data;
    task = task{day}{ep};

    % use the coordinate program to fetch the track coordinates
    coordInTime = task.linearcoord;% 4 segs * 2 (x, y) * time (see below for seg)
    % more about coordInTime: see getcoord_wtrack
    %   1    2    3
    %   |    |    |
    %   |    |    |
    %   |    |    |
    %   |    |    |
    %   4----5----6
    % coordInTime{1}, left traj (aka,[2 5 4 1]); coordInTime{2}, right traj (aka,[2 5 6 3])
    clear task;

    % If sj_velocitydayprocess resulted in couple of extra points in pos
    % structure, and your task file hasnt been updated since, skip first few
    % points (usually 2)

    % match coordInTime and pos time. Only deal with the case that pos time is
    % longer
    if size(pos,1)~=length(coordInTime{1})
        rem = abs(size(pos,1) - length(coordInTime{1}));
        if size(pos,1)>length(coordInTime{1})
            pos=pos(rem+1:end,:);   
        end
    end

    %%
    %initialize variables
    timestep = pos(2,1) - pos(1,1);
    poslength = size(pos,1);
    segment = zeros(poslength,1);
    lindist = ones(size(pos,1), 1) * -1; % -1 as initial value
    vect = zeros(poslength,2);
    newpos = pos;
    %%
    %calculate the linear distance to each coordinate. Distance from the reference, the Center Well.
    distsum = [];% dis: 2to5,2to4, 2to1 as sqrt(x^2 + y^2)
    for i = 1:length(coordInTime)%2 trajs
        firstcoord{i} = coordInTime{i}(:,:,1);
        if (size(firstcoord{i},1) > 1)
            distsum{i}(1,1) = 0;
            for j = 2:size(firstcoord{i},1)% 4 segs
                distsum{i}(j,1) = sqrt( ((firstcoord{i}(j,1) - firstcoord{i}(j-1,1))^2) + ((firstcoord{i}(j,2) - firstcoord{i}(j-1,2))^2) ) + ...
                    distsum{i}(j-1);
            end
        end
    end
    %%
    % getSegmentTable
    wells = [];
    wellindex = [];
    trajwells = [];

    tmpcoord = [];
    coordtraj = [];
    trajlength = []; 

    for i = 1:length(firstcoord)% i = 1
        if ~rowfind(firstcoord{i}(1,:),wells)
            wellindex = [wellindex; [i 1]];% wellindex 1 = coord 2
            wells = [wells; firstcoord{i}(1,:)];% well 1: (x, y) of coord 2
            trajwells(i,1) = size(wellindex,1);
        else
            trajwells(i,1) = rowfind(firstcoord{i}(1,:),wells);
        end
        if ~rowfind(firstcoord{i}(size(firstcoord{i},1),:),wells)%well 2, last coord (aka 1)
           wellindex = [wellindex; [i size(firstcoord{i},1)]];
           wells = [wells; firstcoord{i}(size(firstcoord{i},1),:)];
           trajwells(i,2) = size(wellindex,1);
        else
           trajwells(i,2) = rowfind(firstcoord{i}(size(firstcoord{i},1),:),wells);
        end

        tmpcoord = [tmpcoord ; [firstcoord{i}(1:end-1,:) firstcoord{i}(2:end,:)]]; %get the coordinates of the start and end points of each segment
        % tmpcoord: 3 segs* [(x,y) (x,y)]:2to5, 5to4, 4to1
        tmpcoordtraj = ones(size(firstcoord{i},1)-1,1)*i;
        tmpcoordseg = (1:(size(firstcoord{i},1)-1))';
        trajlength(i,1:2) = [size(firstcoord{i},1)-1 size(coordtraj,1)];
        coordtraj = [coordtraj;[tmpcoordtraj tmpcoordseg]];
    end


    segmentCoords = [];
    coordind = [];
    indcount = 0;
    for i = 1:size(tmpcoord,1)% segment loop
        [test, testind] = ismember(tmpcoord(i,:),segmentCoords,'rows','legacy');
        reverse = 0;
        if ~(test)
            %if the segment wasn't recorded in one direction, check the reverse
            [test, testind] = ismember([tmpcoord(i,3:4) tmpcoord(i,1:2)],segmentCoords,'rows','legacy');
            reverse = 1;
        end

        origin = tmpcoord(i,1:2);% start point
        endpoint = tmpcoord(i,3:4);
        seglength = sqrt( ((endpoint(1) - origin(1))^2) + ((endpoint(2) - origin(2))^2) ); % length
        if ~(test)
            %this segment hasn't been added to the new list yet, so give it a
            %new index
            segmentCoords = [segmentCoords; tmpcoord(i,:)];
            indcount = indcount+1;
            segmentLength(indcount) = seglength;
            coordind = [coordind; indcount];
            %record which trajectories this segment belongs to
            segmentTrajectories(indcount,1) = bitset(0,coordtraj(i,1),1); % Middle arm belongs to traj 1 and 2, so becomes 3                
        else
            %this segment is already in the new list, so give it the original
            %index
            coordind = [coordind;testind];
            segmentTrajectories(testind) = bitset(segmentTrajectories(testind),coordtraj(i,1),1);
        end
    end
    segmentInfo.segmentCoords = segmentCoords;% 5segs *[(x,y) (x,y)]: 2to5 5to4 4to1, 5to6, 6to3
    segmentInfo.segmentLength = segmentLength;% length for 5 segs
    segmentInfo.segmentTraj = segmentTrajectories; % Still oly 2 traj: 1 and 2, 3 is for center arm
    segmenttable = [coordtraj coordind];% 1: left or right; 2: 3 segs; 3: segs in the tabel (1-5)
    %%
    %find which segments connect to the start and end of each segment
    % more about segs:

    %   |    |    |
    %   |    |    |
    %   |3   | 1  |5
    %   |    |    |
    %   -----------
    %     2     4
    for i = 1:size(segmentCoords,1)
        tmp = find( ((segmentCoords(:,1) == segmentCoords(i,1))&(segmentCoords(:,2) == segmentCoords(i,2))) | ((segmentCoords(:,3) == segmentCoords(i,1))&(segmentCoords(:,4) == segmentCoords(i,2))) );
        startLinkSegments{i} = setdiff(tmp,i);   
        tmp = find( ((segmentCoords(:,1) == segmentCoords(i,3))&(segmentCoords(:,2) == segmentCoords(i,4))) | ((segmentCoords(:,3) == segmentCoords(i,3))&(segmentCoords(:,4) == segmentCoords(i,4))) );
        endLinkSegments{i} = setdiff(tmp,i);  
    end
    % 2 links 1,3; 3 links 2; 4 links 1,2; 5 links 4.
    %%
    connectivityTable = [];
    wellsegments = [];
    uniqueSegments = unique(segmenttable(:,3))';
    distanceDivisor = 1000;
    if (length(uniqueSegments) > 1)
        for i = uniqueSegments
            segindex = find(segmenttable(:,3) == i)';
            for j = segindex
                traj = segmenttable(j,1);
                segnum = segmenttable(j,2);
                seglength = segmentLength(i);
                %find the segment indeces of the trajectory endpoints
                if (sum((wellindex(:,1) == traj)&((wellindex(:,2) == segnum)|(wellindex(:,2) == (segnum+1)))))
                    transindex = min(find((wellindex(:,1) == traj)&((wellindex(:,2) == segnum)|(wellindex(:,2) == (segnum+1)))));
                    wellsegments = [wellsegments; [i transindex]];
                end

                %find the segment indeces of all segments that are directly
                %connected to this segment             
                thisSegCoord = segmentCoords(i,:);
                connectedSegments = find( ((segmentCoords(:,1) == thisSegCoord(1))&(segmentCoords(:,2) == thisSegCoord(2))) | ...
                                           ((segmentCoords(:,1) == thisSegCoord(3))&(segmentCoords(:,2) == thisSegCoord(4))) | ...
                                           ((segmentCoords(:,3) == thisSegCoord(1))&(segmentCoords(:,4) == thisSegCoord(2))) | ...
                                           ((segmentCoords(:,3) == thisSegCoord(3))&(segmentCoords(:,4) == thisSegCoord(4))) ); 

                connectedSegments = unique(setdiff(connectedSegments,i));
                %connectedSegments = segmenttable(find((segmenttable(:,1) == traj) & ( (segmenttable(:,2) == segnum+1)|(segmenttable(:,2) == segnum-1) )),3);

                distanceDivisor = 1000;
                for k = connectedSegments'
                    %the (i,k)th entry in the connectivity matrix says if the
                    %ith segment is directly connected to the kth segment,
                    %and by what distance: (exp(-lengthOfSegmentI))
                    connectivityTable(i,k) = exp(-seglength/distanceDivisor); 
                end
            end
        end
    else
        %there is only one segment
        wellsegments = [1 1];
        connectivityTable = exp(-segmenttable(1,3));
    end
    % more about connectivityTable, aka ct: 5segs * 5 segs, if seg i links seg j:
    % ct(i,j) = seglength(i)
    wellsegments = sortrows(wellsegments,2);
    wellsegments = wellsegments(:,1);% 1,3,5 segs link wells


    for i = 1:length(wellsegments)
        %calculate the length of the arms leading up to the wells (until an
        %intersection is hit). This will be used later to determine if the
        %animal completed a trajectory
        armlength(i) = 0;
        segcount = 1;
        foundIntersection = 0;
        segindex = wellsegments(i);
        tempTable = connectivityTable;
        while( (foundIntersection == 0) & (segcount < length(uniqueSegments)) )
            if(sum(tempTable(segindex,:) > 0) > 1)
                foundIntersection = 1;            
            else
                %we have not reached an intersection yet, so find the next
                %segments in the connection tree
                tempTable = tempTable*connectivityTable;
                segcount = segcount + 1;
            end
        end
        %calculate the length of the independant arm (all nonzero numbers on
        %the row should be the same, so just pick the maximum one)           
        armlength(i) = -log(max(tempTable(segindex,:)));   % interesting calculation. use e^(x+y) = e^x*e^y   
    end
    % more about armlength: 3 well segs: 1, 3, 5. length = 1;2+3;4+5

    %%
    tempTable = connectivityTable;
    distanceTable = connectivityTable;
    pathTable = cell(size(connectivityTable));
    diagIndeces = find(eye(size(connectivityTable,1))); %these are the indeces to the diagonal entries
    [nonZerosi,nonZerosj] = find(distanceTable > 0);
    for i = 1:length(nonZerosi)
        pathTable{nonZerosi(i),nonZerosj(i)} = [nonZerosi(i) nonZerosj(i)];
    end

    while(sum(distanceTable(:) == 0) > 0)   %we do this loop until we have found a path between every pair of segments 


        oldtempTable = tempTable;
        tempTable = tempTable*connectivityTable;  %this finds the grandchildren...greatgrandchildren...and so on, of each segment, and the corresponding distance 
        zerovals = find(distanceTable == 0);
        %which segment pairs became connected?
        [switchedValsi, switchedValsj] = find( (distanceTable == 0) & (tempTable > 0) );

        %we only update the values that were zero (because we already found the
        %minimum path for the other pairs)
        distanceTable(zerovals) = tempTable(zerovals);

        for i = 1:length(switchedValsi)
            if (switchedValsi(i) ~= switchedValsj(i))
                %find which segment links one path to the new segment

                linker = find( (oldtempTable(switchedValsi(i),:) > 0) & (connectivityTable(:,switchedValsj(i))' > 0) );
                prepath = pathTable{switchedValsi(i),linker}(1:end-1);
                postpath = pathTable{linker,switchedValsj(i)}(2:end);
                %add the path to pathTable
                pathTable{switchedValsi(i),switchedValsj(i)} = [prepath linker postpath];
            end
        end 
    end
    for i = 1:length(diagIndeces)
        %the diagonal entries are wrong, so we fix them
        distanceTable(diagIndeces(i)) = 1;
        pathTable{diagIndeces(i)} = i;
    end

    %convert the exponent distance back to normal distance
    %distanceTable =(-log(distanceTable)) - log(multtable);
    distanceTable =-log(distanceTable) * distanceDivisor; %path distance

    %only keep the distances from the wells (we don't care about the other
    %segment pairs)
    distanceTable = distanceTable(wellsegments,:);
    pathTable = pathTable(wellsegments,:);
    segmentDirection = [];
    % more about distanceTable: path 1-2-3, distance = 1+2; 
    %%
    %finally, we need to calculate whether a segment is aligned in the foreward
    %or backward direction relative to each well
    for i = 1:size(pathTable,1)
        for j = 1:size(pathTable,2)        
            if ( (isempty(startLinkSegments{j})) & (isempty(endLinkSegments{j})) )
                %there is only one segment
                if (i==1)
                    segmentDirection(i,j) = 1;
                else
                    segmentDirection(i,j) = 0;
                end              
            else
                foundpath = 0;
                for k = pathTable{i,j}
                    if (ismember(k,startLinkSegments{j}))
                        segmentDirection(i,j) = 1; %from this well, the segment will linearize in the foreward direction
                        foundpath = 1;
                        break;
                    elseif (ismember(k,endLinkSegments{j}))
                        segmentDirection(i,j) = 0; %or from the backward direction
                        foundpath = 1;
                        break;
                    end
                end
                if (foundpath == 0)
                    %this is the same segment as the well, and it will
                    %therefore have either no start segments or no end
                    %segments
                    if (isempty(startLinkSegments{j}))
                        segmentDirection(i,j) = 1;% forward: center to outer, 1
                    elseif (isempty(endLinkSegments{j}))
                        segmentDirection(i,j) = 0;% backward: outer to center, 0
                    end
                end            
            end
        end
    end

    %create the wellSegmentInfo structure          
    wellSegmentInfo.distanceTable = distanceTable;
    wellSegmentInfo.segmentIndex = wellsegments;
    wellSegmentInfo.distanceToIntersection = armlength; 
    wellSegmentInfo.wellCoord = wells;
    wellSegmentInfo.pathTable = pathTable;
    wellSegmentInfo.segmentDirection = segmentDirection; 
%}
%% Potential and Gradient
%%%%%%%%%%%%%%%%%%%%%%%%%

seg_start = segmentCoords(:,1:2);
seg_end = segmentCoords(:,3:4);
seg_length = segmentLength;

segmentCoords = round(segmentCoords,1);
seg_dirs = segmentCoords(:,3:4) - segmentCoords(:,1:2);
for i=1:size(seg_dirs,1)
    seg_dirs(i,:) = seg_dirs(i,:)/norm(seg_dirs(i,:));
    
end

%%%%%%%%%%%%%%%%%
%Backup variables (will be used at end)
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
    seg_dirs = seg_dirs([1 2 3 5],:);
    seg_start = seg_start([1 2 3 5],:);
    seg_length = seg_length([1 2 3 5]);
    %segmentCoords = segmentCoords([1 2 3 5],:);
end
%}


%%

raw_pos= round([pos(:,2) pos(:,3)],1);
buffer = 20; % how much extra around edges after max & min pos.
den = 0.1;

[X,Y] = meshgrid((min(raw_pos(:,1))-buffer):den:(max(raw_pos(:,1))+buffer),(min(raw_pos(:,2))-buffer):den:(max(raw_pos(:,2))+buffer));  % covers all points with grid spacing = den
X = round(X,1);
Y = round(Y,1);
pots = zeros(size(X,1),size(X,2),5); % for 5 segments
w1 = 1/50;  %Wider contours avoid sharp changes in projection direction for points at the intersections, 1/50 works well
w2 = 1/1e6;

points = [];
for i = [1 2 3 5] %1:size(segmentCoords,1)
    if ~mod(i,2) % even segment = horizontal
        coeffs = polyfit(segmentCoords(i,[1 3]),segmentCoords(i,[2 4]),1); % fit for y = mx+c
        pots(:,:,i) = (((1+exp(-1.*((((coeffs(1).*X)+coeffs(2))-Y).^2).*w1)).^-1)-1)+...
                       (((1+exp(-1.*((((coeffs(1).*X)+coeffs(2))-Y).^2).*w2)).^-1)-1) ;
        %pots(
    else         % odd segment = vertical
        coeffs = polyfit(segmentCoords(i,[2 4]),segmentCoords(i,[1 3]),1);  % fit for x = my+c
        pots(:,:,i) = (((1+exp(-1.*((((coeffs(1).*Y)+coeffs(2))-X).^2).*w1)).^-1)-1)+...
                        (((1+exp(-1.*((((coeffs(1).*Y)+coeffs(2))-X).^2).*w2)).^-1)-1);       
        points(i,:) = segmentCoords(i,[1,3]);
    end
end
points = sort(points,2);
f_pots = sum(pots,3);


%points = [77.2 13.5; 31.6 10.8; 120.8 16.1];
f = @(x,y,x1,y1,w) (1+exp(-(w1.*(x-x1).^2+w1.*(y-y1).^2))).^-1;

%f_pots2 = f_pots;
%
tic;
for i = [1 3 5]
    min_points = find(f_pots == min(f_pots(X>(points(i,1)-10) & X<(points(i,2)+10))));
    x_linspace = X(min_points);
    y_linspace = Y(min_points);
    
    for j = 1:length(x_linspace)
    z =feval(f,X,Y,x_linspace(j),y_linspace(j),w1);
    f_pots = f_pots - z;
    end
disp(i)
end

%}

[px,py] = gradient(f_pots); % Calculate gradient

%{
figure
contour(X,Y,f_pots)
s = surf(X,Y,f_pots)
set(s,'LineStyle','none')
hold on

quiver(X,Y,px,py);
scatter(raw_pos(:,1),raw_pos(:,2),5)
hold off
%}


%% Projection

%% 

%Calculate how X and Y change with index of matrix since meshgrid creates
%linear spacing. Use indices to lookup values for gradient for every
%position point.

x_fit = polyfit(X(1,1:10),1:10,1);
y_fit = polyfit(Y(1:10,1)',1:10,1);

f_grad = [];
for i=1:length(raw_pos(:,1))
    %x index gives columns and y index gives rows
    x_ind = int16(x_fit(1)*raw_pos(i,1) + x_fit(2)); 
    y_ind = int16(y_fit(1)*raw_pos(i,2) + y_fit(2));
    
    f_grad(i,:) = [px(y_ind,x_ind) py(y_ind,x_ind)];
    f_grad(i,:) = f_grad(i,:)/norm(f_grad(i,:));
    
end


%%
%Calculate projection of point along respective gradient by changing basis
%and using unit vector along gradient.

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
    
    %Select closest segment
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


%%
%%%%%%%%%%%%%%%%%%%%
%lindist calculation
%%%%%%%%%%%%%%%%%%%%
%1)Index of segment projection is on new_pos
%2)Is segment start point for  wards or backwards of all of the 3 wells
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

%Adjust connectivity table from old code.
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