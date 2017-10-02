load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, 10}.trajtime;

names = cell(1,size(traj_time,1));
names(:) = {'fig'};
file_ext = cell(1,size(traj_time,1));
file_ext(:) = {'.png'};
names = matlab.lang.makeUniqueStrings(names);
names = strcat(names,file_ext);


for i = 1:10 %%size(traj_time,1)
    
    
    
    for  lr = 1:2
        start_time = trajtimes{lr}(i,1);
        end_time = trajtimes{lr}(i,2);
        start_ind = find(traj_pos{lr}(:,3)==start_time);
        end_ind = find(traj_pos{lr}(:,3)==end_time);
    
        figure
        contour(X,Y,f_pots(:,:,lr))
        hold on
        scatter(traj_pos{lr}(start_ind:end_ind,1),traj_pos{lr}(start_ind:end_ind,2))
        scatter(new_pos{lr}(start_ind:end_ind,1),new_pos{lr}(start_ind:end_ind,2))
        plot(new_pos{lr}(start_ind:end_ind,1),new_pos{lr}(start_ind:end_ind,2))
        q = quiver(traj_pos{lr}(start_ind:end_ind,1),traj_pos{lr}(start_ind:end_ind,2),-quiv{lr}(start_ind:end_ind,1),-quiv{lr}(start_ind:end_ind,2));
        q.AutoScaleFactor = 0.5;
        pause(0.0001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1); 
    end
end