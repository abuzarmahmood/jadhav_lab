load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, ep}.trajtime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For traj divided projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:10 %%size(traj_time,1)
    
    
    
    for  lr = 1:2
        start_time = traj_time{lr}(i,1);
        end_time = traj_time{lr}(i,2);
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
%%%%%%%%%%%%%%%%%%%%%%
%For potential field 2
%%%%%%%%%%%%%%%%%%%%%%
load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, ep}.trajtime;
wells_tend = trajinfo{1, 1}{1, ep}.wellstend;

for i = 1:size(traj_time,1)

        start_time = traj_time(i,1);
        end_time = traj_time(i,2);
        start_ind = find(pos(:,1)==start_time);
        end_ind = find(pos(:,1)==end_time);
    
        figure
        contour(X,Y,f_pots)
        hold on
        scatter(pos(start_ind:end_ind,2),pos(start_ind:end_ind,3))
        scatter(new_pos(start_ind:end_ind,1),new_pos(start_ind:end_ind,2))
        plot(new_pos(start_ind:end_ind,1),new_pos(start_ind:end_ind,2))
        q = quiver(pos(start_ind:end_ind,2),pos(start_ind:end_ind,3),-f_grad(start_ind:end_ind,1),-f_grad(start_ind:end_ind,2));
        q.AutoScaleFactor = 0.5;
        pause(0.0001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1); 
        pbaspect([1 1 1])
    end
end


%%
%%%%%%%%%%%%%%%%
%For Wenbo linpos
%%%%%%%%%%%%%%%%%
close all;
load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, ep}.trajtime;
wells_tend = trajinfo{1, 1}{1, ep}.wellstend;

for i = 1:size(traj_time,1)
    
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
        plot(pos(start_ind:end_ind,1),lindist(start_ind:end_ind,wells_tend(i,1)))
        
        %plot(pos(start_ind:end_ind,1),lindist(start_ind:end_ind,1:3))
end

figure
scatter(pos(:,2),pos(:,3))
hold on
plot(newpos(:,2),newpos(:,3))
hold off

%%
%%%%%%%%%%%%%%%%
%For Abu linpos
%%%%%%%%%%%%%%%%%
close all;
load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, ep}.trajtime;
wells_tend = trajinfo{1, 1}{1, ep}.wellstend;

for i = 1:size(traj_time,1)
    
        start_time = traj_time(i,1);
        end_time = traj_time(i,2);
        start_ind = find(pos(:,1)==start_time);
        end_ind = find(pos(:,1)==end_time);
    
        figure
        subplot(2,1,1)
        hold on
        scatter(raw_pos(start_ind:end_ind,1),raw_pos(start_ind:end_ind,2))
        scatter(new_pos(start_ind:end_ind,1),new_pos(start_ind:end_ind,2))
        plot(new_pos(start_ind:end_ind,1),new_pos(start_ind:end_ind,2))
        
        subplot(2,1,2)
        hold on
        plot(pos(start_ind:end_ind,1),lindist_abu(start_ind:end_ind,wells_tend(i,1)))
        
        %plot(pos(start_ind:end_ind,1),lindist(start_ind:end_ind,1:3))
end

figure
scatter(pos(:,2),pos(:,3))
hold on
plot(new_pos(:,1),new_pos(:,2))
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Side-by-side comparison (potential field 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;
load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, ep}.trajtime;
%traj_time = [2180 2195;2070 2085; 2095 2100; 2245 2270; 2395 2410; 2670 2700; 2252.8 2253];
%traj_time = [2252.8 2253];
wells_tend = trajinfo{1, 1}{1, ep}.wellstend;

for i = 1:size(traj_time,1)

        %start_time = traj_time(i,1);
        %end_time = traj_time(i,2);
        %start_ind = find(pos(:,1)==start_time);
        %end_ind = find(pos(:,1)==end_time);
        [start_val start_ind] = min(abs(pos(:,1)-traj_time(i,1)));
        [end_val end_ind] = min(abs(pos(:,1)-traj_time(i,2)));
        
        l = end_ind-start_ind+1;
        
        c = linspace(1,10,l);
        
        
        figure
        subplot(2,2,1)
        hold on
        scatter(pos(start_ind:end_ind,2),pos(start_ind:end_ind,3))
        scatter(newpos(start_ind:end_ind,2)+linspace(0,5,l)',newpos(start_ind:end_ind,3)+linspace(0,5,l)',[],c)
        plot(newpos(start_ind:end_ind,2),newpos(start_ind:end_ind,3))
        colormap(hsv)
        
       
        subplot(2,2,2)
        hold on
        scatter(pos(start_ind:end_ind,2),pos(start_ind:end_ind,3))
        %scatter(new_pos(start_ind:end_ind,1)+linspace(0,5,l)',new_pos(start_ind:end_ind,2)+linspace(0,5,l)',[],c)
        scatter(new_pos(start_ind:end_ind,1),new_pos(start_ind:end_ind,2),[],c)
        plot(new_pos(start_ind:end_ind,1),new_pos(start_ind:end_ind,2))
        %contour(X,Y,f_pots)
        colormap(hsv)
        %q = quiver(pos(start_ind:end_ind,2),pos(start_ind:end_ind,3),-f_grad(start_ind:end_ind,1),-f_grad(start_ind:end_ind,2));
        %q.AutoScaleFactor = 0.5;
        pause(0.0001);
        
        
        subplot(2,2,3)
        %plot(pos(start_ind:end_ind,1),lindist(start_ind:end_ind,wells_tend(i,1)))
        plot(start_ind:end_ind,lindist(start_ind:end_ind,wells_tend(i,1)))
        
        subplot(2,2,4)
        plot(pos(start_ind:end_ind,1),lindist_abu(start_ind:end_ind,wells_tend(i,1)))
        
        %frame_h = get(handle(gcf),'JavaFrame');
        %set(frame_h,'Maximized',1); 
        
end

figure
plot(newpos(:,2),newpos(:,3))
hold on
plot(new_pos(:,1),new_pos(:,2))
hold off


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distance from raw position
%%%%%%%%%%%%%%%%%%%%%%%%%%%
orig_difference = pos(:,2:3)-newpos(:,2:3);
abu_difference = pos(:,2:3)-new_pos(:,1:2);

for i = 1:size(pos,1)
    orig_lin_diff(i) = norm(orig_difference(i,:));
    abu_lin_diff(i) = norm(abu_difference(i,:));
end

figure
plot(orig_lin_diff)
hold on
plot (abu_lin_diff)
dim = [.2 .5 .3 .3];
txt = ['Original Mean %0.2f \t SD %0.2f'  newline  'Modified Mean %0.2f \t SD %0.2f'];
str = sprintf(txt,mean(orig_lin_diff(orig_lin_diff<20)),std(orig_lin_diff(orig_lin_diff<20)),mean(abu_lin_diff),std(abu_lin_diff));
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%%
%%%%%%%%%%%%%
%Linpos movie
%%%%%%%%%%%%%

close all;
load('KL8trajinfo01.mat')
traj_time = trajinfo{1, 1}{1, ep}.trajtime;
wells_tend = trajinfo{1, 1}{1, ep}.wellstend;

v = VideoWriter('linpos_vid.avi');
open(v);

for i = 1:size(traj_time,1)
        figure
        start_time = traj_time(i,1);
        end_time = traj_time(i,2);
        start_ind = find(pos(:,1)==start_time);
        end_ind = find(pos(:,1)==end_time);
    
        nframe=end_ind-start_ind;
        traj_inds = start_ind:end_ind;
        mov(1:nframe)= struct('cdata',[],'colormap',[]);
        set(gca,'nextplot','replacechildren')
        for k=1:nframe
                subplot(2,1,1)
                hold on
                scatter(raw_pos(traj_inds(k),1),raw_pos(traj_inds(k),2),2,[1 0 0],'filled')
                %scatter(pos(traj_inds(k),2),pos(traj_inds(k),3),2,[1 0 0],'filled')
                scatter(new_pos(traj_inds(k),1),new_pos(traj_inds(k),2),2,[0 0 1],'filled')
                %quiver(raw_pos(traj_inds(k),1),raw_pos(traj_inds(k),2),f_grad(traj_inds(k),1),f_grad(traj_inds(k),2))
                %scatter(newpos(traj_inds(k),2),newpos(traj_inds(k),3),2,[0 0 1],'filled')
                axis([min(raw_pos(:,1)) max(raw_pos(:,1)) min(raw_pos(:,2)) max(raw_pos(:,2))])
                %axis([min(pos(:,2)) max(pos(:,2)) min(pos(:,3)) max(pos(:,3))])
                subplot(2,1,2)
                plot(pos(traj_inds(1:k),1),lindist_abu(traj_inds(1:k),wells_tend(i,1)))
                %plot(pos(traj_inds(1:k),1),lindist(traj_inds(1:k),wells_tend(i,1)))
                hold off
              %mov(k)=getframe(gcf);
              writeVideo(v,getframe(gcf))
              
        end
        
        
        
        %plot(pos(start_ind:end_ind,1),lindist(start_ind:end_ind,1:3))
end
close all
close(v)