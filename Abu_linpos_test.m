load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, 10}.trajtime;

names = cell(1,size(traj_time,1));
names(:) = {'fig'};
file_ext = cell(1,size(traj_time,1));
file_ext(:) = {'.png'};
names = matlab.lang.makeUniqueStrings(names);
names = strcat(names,file_ext);


for i = 1:10 %%size(traj_time,1)
    start_time = traj_time(i,1);
    end_time = traj_time(i,2);
    start_ind = find(pos(:,1)==start_time);
    end_ind = find(pos(:,1)==end_time);
    
    figure
    hold on
    contour(X,Y,f_pots)
    scatter(raw_pos(start_ind:end_ind,1),raw_pos(start_ind:end_ind,2))
    scatter(new_pos(start_ind:end_ind,1),new_pos(start_ind:end_ind,2))
    q = quiver(raw_pos(start_ind:end_ind,1),raw_pos(start_ind:end_ind,2),-f_grad(start_ind:end_ind,1),-f_grad(start_ind:end_ind,2));
    q.AutoScaleFactor = 0.5;
    pause(1);
end