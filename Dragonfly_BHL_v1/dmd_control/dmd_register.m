function dmd_register(app)


x_center =601;
y_center = 352;
ang_offset = -.08;
ang_i = ((2/4:1/1:2)+ang_offset)*pi;
ang_o = ((3/2:1/3:2)+ang_offset)*pi;
x0 = [round(cos(ang_o')*200); round(cos(ang_i')*100)]+x_center;
y0 = [round(sin(ang_o')*200); round(sin(ang_i')*100)]+y_center;
% x0 = [
%     761
%     564
%     485
%     697
%     385
%     358
%     ];
% y0 = [
%     440
%     584
%     146
%     356
%     504
%     265
%     ];

% x_o=200;y_o=300;

selectedButton = app.DMDControlButtonGroup.SelectedObject;
switch selectedButton.Text
    case 'Project Pattern'

%         regist_pat = register_pattern_test(app.dmd.device.height,app.dmd.device.width,[x_o,y_o]);
        regist_pat = false(app.dmd.device.height,app.dmd.device.width);
        regist_pat(sub2ind(size(regist_pat),x0,y0))=true;
        app.dmd.project(double(regist_pat))

    case 'Register'
        labels = strsplit(sprintf('%g,',1:length(x0)),','); labels = labels(1:end-1);
        figure(758);plot(x0,y0,'o');text(x0+10,y0+10,labels,'fontsize',12);axis equal
        drawnow;pause(0.1);
        f = dmd_load_reference_img(app,'register');
%         a = f.Children;
%         c = a.Children;
%         im = double(c(2).CData);
%         c.CData = imgaussfilt(im,1);
        colormap(jet);%set(gca,'clim',[prctile(im,10,'all') prctile(im,99.995,'all')])
%         title('Please click outer corner of "L" in the sequence of 00 -> 10 -> 01. Press Enter to finish')
        [xi,yi] = getpts;
%         [xi, yi] = (getline(gca));
        if length(xi)<3,return;end
        
        Xs = [y0(2:end)-y0(1) x0(2:end)-x0(1)]';
        Xc = [xi(2:end)-xi(1) yi(2:end)-yi(1)]';
        reg_rot_dil = Xc/Xs; 
        reg_trans = [xi(1) yi(1)]'-reg_rot_dil*[y0(1) x0(1)]';

%         xd00 = [x_o,y_o]'; xd01 = [x_o app.dmd.device.height-y_o]'; xd10 = [app.dmd.device.width-x_o y_o]';
%         xc00 = [xi(1),yi(1)]'; xc10 = [xi(2),yi(2)]'; xc01 = [xi(3),yi(3)]';
%         Xd_prm = [xd00-xd10 xd10-xd01]; Xc_prm = [xc00-xc10 xc10-xc01];
%         reg_rot_dil = Xc_prm/Xd_prm; reg_trans = xc00-reg_rot_dil*xd00;
%         vec = round([xi,yi]);
        app.registration_data = [reg_trans reg_rot_dil]';
        dlmwrite('dmd_registration_data.txt',app.registration_data,'precision',16);
        disp 'Registration Success!'
end


end