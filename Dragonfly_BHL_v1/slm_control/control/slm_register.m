function slm_register(app)

% x0 = [640.250000000000;1415.75000000000;1438.25000000000]; % full screen
% y0 = [397.250000000000;385.250000000000;739.250000000000]; % full screen
% x0 = [640.250000000000;1415.75000000000;1438.25000000000]; 
% y0 = [500; 490; 600];

% x0 = [497.195652173913;1421.41925465838;491.232919254658]; % old
% y0 = [512.698757763975;534.164596273292;641.493788819876];
% x_center = 1.2164e+03;
% y_center = 877.5815;

% ---------- 125 mm lens --------
x0 = [814.482165605096;783.664331210191;864.377707006369;923.078343949044]+370;
y0 = [792.224840764331;853.860509554140;825.977707006369;858.263057324841]+45;
% -------------------------------

% ---------- 75 mm lens --------
% x0 = [814.482165605096;783.664331210191;864.377707006369;923.078343949044]+170;
% y0 = [792.224840764331;853.860509554140;825.977707006369;858.263057324841]+150;
% -------------------------------

% ang = (2/4:1/4:(2-1/4))*pi;
% x0 = [round(cos(ang')*100)]+x_center;
% y0 = [round(sin(ang')*50)]+y_center;
z0 = zeros(size(x0));

xyz_origin = mean([x0 y0 z0],1);

app.xyz_origin = xyz_origin;

xyz_full = apply_optimal_offsets([x0 y0 z0]',xyz_origin','full');



selectedButton = app.SLMControlButtonGroup.SelectedObject;
switch selectedButton.Text
    case 'Project Pattern'
%         try
%             load('slm_control\control\regist_pat.mat','phase')
%         catch
%             regist_pat = gs_target_gen([x0 y0],'spt');
            phase = wgs_spots(xyz_full',30,1);
%             phase = gs(regist_pat,200);
%             save('slm_control\control\regist_pat.mat','phase')
%         end
        app.SLM.project(phase)
    case 'Register'
        figure;plot(x0,y0,'o');text(x0+10,y0+10,cellstr(sprintf('%g',1:length(x0))'),'fontsize',12)
        fig = slm_load_reference_img(app,'register'); if fig==false;return;end
        [xi,yi] = getpts;
%         [xi, yi] = (getline(gca));
        if length(xi)<3,return;end

%         xd00 = [x0(1) y0(1)]'; xd10 = [x0(2) y0(2)]'; xd01 = [x0(3) y0(3)]';
%         xc00 = [xi(1),yi(1)]'; xc10 = [xi(2),yi(2)]'; xc01 = [xi(3),yi(3)]';
%         Xd_prm = [xd00-xd10 xd10-xd01]; Xc_prm = [xc00-xc10 xc10-xc01];
%         reg_rot_dil = Xc_prm/Xd_prm; reg_trans = xc00-reg_rot_dil*xd00;
%         vec = round([xi,yi]);

        xyz = apply_optimal_offsets([x0 y0 z0]',xyz_origin','apply_mat');
        Xs = xyz(:,2:end)-xyz(:,1);
        Xc = [xi(2:end)-xi(1) yi(2:end)-yi(1) ones(length(x0)-1,1)]';
        reg_rot_dil = Xc/Xs; reg_trans = [xi(1) yi(1) 1]'-reg_rot_dil*xyz(:,1);

        app.registration_data = [reg_trans reg_rot_dil]';
        dlmwrite('slm_registration_data.txt',app.registration_data,'precision',20);

%          xyz = apply_optimal_offsets([x0 y0 z0]',xyz_origin','apply_mat');
%         Xs = [xyz(:,2)-xyz(:,1) xyz(:,3)-xyz(:,1) xyz(:,4)-xyz(:,1)];
%         Xc = [xi(2)-xi(1) yi(2)-yi(1) 1; xi(3)-xi(1) yi(3)-yi(1) 1; xi(4)-xi(1) yi(4)-yi(1) 1]';
%         reg_rot_dil = Xc/Xs; reg_trans = [xi(1) yi(1) 1]'-reg_rot_dil*xyz(:,1);
%         
%         app.registration_data = [reg_trans reg_rot_dil]';
%         dlmwrite('slm_registration_data.txt',app.registration_data);
        disp 'Registration Success!'
    case 'Live Register'
        df_app = evalin('base','hDragonflyApp');
        [cam_r,cam_c] = ind2sub(size(df_app.camApp.cameraMaskImg),find(df_app.camApp.cameraMaskImg));
        cam_offset = [cam_c(1) cam_r(1)]-1;
%         try
%             load('slm_control\control\live_regist_pat.mat','phase_all')
%         catch
            xi = zeros(size(x0));
            yi = zeros(size(y0));
            for ii = 1:length(x0)
%                 regist_pat = gs_target_gen([x0(ii) y0(ii)],'spt');
%                 phase_all(:,:,ii) = gs(regist_pat,200);
                phase_all(:,:,ii) = wgs_spots(xyz_full(:,ii)',30,1);
            end
%             save('slm_control\control\live_regist_pat.mat','phase_all')
%         end
        
        for ii = 1:length(x0)
            app.SLM.project(squeeze(phase_all(:,:,ii)))
            pause(.5)
            live_im=df_app.camApp.snapOneFrame;
            figure(998);clf;imshow(live_im,[],'initialmagnification','fit');colormap(jet)
            [xii,yii] = getpts;

            xi(ii) = xii+cam_offset(1);
            yi(ii) = yii+cam_offset(2);
        end
        
%         xd00 = [x0(1) y0(1)]'; xd10 = [x0(2) y0(2)]'; xd01 = [x0(3) y0(3)]';
%         xc00 = [xi(1),yi(1)]'; xc10 = [xi(2),yi(2)]'; xc01 = [xi(3),yi(3)]';
%         Xd_prm = [xd00-xd10 xd10-xd01]; Xc_prm = [xc00-xc10 xc10-xc01];
%         reg_rot_dil = Xc_prm/Xd_prm; reg_trans = xc00-reg_rot_dil*xd00;
%         vec = round([xi,yi]);
%         app.registration_data = [reg_trans reg_rot_dil]';
        xyz = apply_optimal_offsets([x0 y0 z0]',xyz_origin','apply_mat');
        Xs = xyz(:,2:end)-xyz(:,1);%[xyz(:,2)-xyz(:,1) xyz(:,3)-xyz(:,1) xyz(:,4)-xyz(:,1)];
        Xc = [xi(2:end)-xi(1) yi(2:end)-yi(1) ones(length(x0)-1,1)]';%[xi(2)-xi(1) yi(2)-yi(1) 1; xi(3)-xi(1) yi(3)-yi(1) 1; xi(4)-xi(1) yi(4)-yi(1) 1]';
        reg_rot_dil = Xc/Xs; reg_trans = [xi(1) yi(1) 1]'-reg_rot_dil*xyz(:,1);

        app.registration_data = [reg_trans reg_rot_dil]';
        dlmwrite('slm_registration_data.txt',app.registration_data,'precision',20);
        
        disp 'Registration Success!'
end


end