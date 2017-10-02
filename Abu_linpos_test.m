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

%%
%%%%%%%%%%%%%%%%
%For Wenbo linpos
%%%%%%%%%%%%%%%%%
close all;
load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, 10}.trajtime;
wells_tend = trajinfo{1, 1}{1, 10}.wellstend;

for i = 1:10 %%size(traj_time,1)
    
        start_time = traj_time(i,1);
        end_time = traj_time(i,2);
        start_ind = find(pos(:,1)==start_time);
        end_ind = find(pos(:,1)==end_time);
    
        figure
        subplot(2,1,1)
        hold on
        scatter(pos(start_ind:end_ind,2),pos(start_ind:end_ind,3))
        scatter(newpos(start_ind:end_ind,2),newpos(start_ind:end_ind,3))
        plot(newpos(start_ind:end_ind,2),newpos(start_ind:end_ind,3))
        
        subplot(2,1,2)
        hold on
        plot(lindist(start_ind:end_ind,wells_tend(i,1)))
end

figure
scatter(pos(:,2),pos(:,3))
hold on
plot(newpos(:,2),newpos(:,3))
hold off