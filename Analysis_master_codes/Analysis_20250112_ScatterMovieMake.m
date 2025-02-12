%continue from
%/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes/Analysis_20241231_Statistics_Spontaneous.m
%line # 811

t_start=-200; t_interval=3;
fig=figure;
set(gcf, 'Position', [100, 100, 1000, 800]); % Set figure window size [left, bottom, width, height]
hold on;
grid on;

   set(gca, 'FontSize', 40); % Set font size for axes

% Prepare video writer
videoWriter = VideoWriter(fullfile(fpath{f},'BasalApical_CSSS_spike'), 'MPEG-4'); % Save as MP4 file
videoWriter.FrameRate = 10; % Adjust frame rate if needed
open(videoWriter);
[V D]=get_eigvector(BASub_silent(:,sum(isnan(BASub_silent))==0),2);
for t=1:400

    [~, BApostSS]=get_STA([Basaltr; Apicaltr], allSpikeClassMat{f}(1,:).*(BlueStim{f}==0),-(t+t_start),t+t_start+t_interval); %pre simple spike
    BApostSS=mean(BApostSS,3,'omitnan');

    [~, BApostCS1st]=get_STA([Basaltr; Apicaltr], allSpikeClassMat{f}(2,:).*(BlueStim{f}==0),-(t+t_start),t+t_start+t_interval);
    BApostCS1st=mean(BApostCS1st,3,'omitnan');

    scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-2,3,100),linspace(-2,3,100)); hold all
    plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);
    l(1)=scatter(BApostCS1st(1,:),BApostCS1st(2,:),80,'filled','Marker','>','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.9); hold all
    l(2)=scatter(BApostSS(1,:),BApostSS(2,:),80,'filled','Marker','o','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.9);
    quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,1)*sqrt(D(1)), V(2,1)*sqrt(D(1)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 1');
    quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,2)*sqrt(D(2)), V(2,2)*sqrt(D(2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 2');
    xlim([-1 2.8])
    ylim([-1.5 2.3])
    colormap(turbo); xlabel('Basal'); ylabel('Apical'); legend(l,{'1st CS','SS'})
    title(['From ' num2str((t+t_start)) ' to ' num2str(t+t_start+t_interval) ' ms'])

    frame = getframe(gcf);
    writeVideo(videoWriter, frame);
    
    cla;
end

% Finalize
hold off;
close(videoWriter);
close(fig);

disp(['Animation saved as ', fullfile(fpath{f},'BasalApical_CSSS_spike')]);