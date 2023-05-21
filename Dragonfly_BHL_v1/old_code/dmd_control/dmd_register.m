function dmd_register(app)

x_o=200;y_o=300;

selectedButton = app.DMDControlButtonGroup.SelectedObject;
switch selectedButton.Text
    case 'Project Pattern'

        regist_pat = register_pattern_test(app.dmd.device.height,app.dmd.device.width,[x_o,y_o]);
        app.dmd.project(double(regist_pat))

    case 'Register'
        f = dmd_load_reference_img(app,'register');
        a = f.Children;
        c = a.Children;
        im = c.CData;
%         colormap(jet);%set(gca,'clim',[prctile(im,10,'all') prctile(im,99,'all')])
        title('Please click outer corner of "L" in the sequence of 00 -> 10 -> 01. Press Enter to finish')
        [xi,yi] = getpts;
%         [xi, yi] = (getline(gca));
        if length(xi)<3,return;end

        xd00 = [x_o,y_o]'; xd01 = [x_o app.dmd.device.height-y_o]'; xd10 = [app.dmd.device.width-x_o y_o]';
        xc00 = [xi(1),yi(1)]'; xc10 = [xi(2),yi(2)]'; xc01 = [xi(3),yi(3)]';
        Xd_prm = [xd00-xd10 xd10-xd01]; Xc_prm = [xc00-xc10 xc10-xc01];
        reg_rot_dil = Xc_prm/Xd_prm; reg_trans = xc00-reg_rot_dil*xd00;
%         vec = round([xi,yi]);
        app.registration_data = [reg_trans reg_rot_dil]';
        dlmwrite('dmd_registration_data.txt',app.registration_data);
        disp 'Registration Success!'
end


end