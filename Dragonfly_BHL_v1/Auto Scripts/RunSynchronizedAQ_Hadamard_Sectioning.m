figX = figure('name','X Motor');
figY = figure('name','Y Motor');

XmotorSerialNumber = 80855716; % x linear actuator
YmotorSerialNumber = 80855718; % y linear actuator

motorX = APTmotor   (XmotorSerialNumber,[0, 0, 500, 300],figX);
motorY = APTmotor   (YmotorSerialNumber,[0, 0, 500, 300],figY);
motorZ = hDragonflyApp.motorZ;

%%
Z0 = motorZ.position; %4.0754
X0 = 18.3298;
Y0 = 10.4843;

Z1 = Z0+.05;
X1 = 16.5870;
Y1 = 12.1913;

zScale = 10/50; % n rpts at z = 50um
nZ = 50;
nX = 3;
nY = 3;

zStep = (Z1-Z0)/nZ;
xStep = 1*sign(X1-X0);
yStep = 1*sign(Y1-Y0);
%%
file_name_prefix = 'M-YQ0201-12';
xSign = true;
ySign = true;
zSign = true;
tic
for ix = 1:nX
%     if xSign
    xi = X0+xStep*(ix-1);
%     else
%     xi = X0+xStep*(nX-ix);        
%     end
%     tic
    motorX.moveTo(xi);

    for iy = 1:nY
        if ySign
            yi = Y0+yStep*(iy-1);
        else
            yi = Y0+yStep*(nY-iy);
        end
%         tic
        motorY.moveTo(yi)
        
        for iz = 1:nZ
%             if zSign
                zi = Z0+zStep*(iz-1);
%             else
%                 zi = Z0+zStep*(nZ-iz);
%             end
%             tic
            motorZ.moveTo(zi)
            file_name = sprintf([file_name_prefix 'x%02u_y%02u_z%02u'],ix,iy,iz);
            hDragonflyApp.appendedfoldernameEditField.Value = file_name;

            nrpt = ceil(iz*zScale);
            hDragonflyApp.RepeatEditField.Value = nrpt;
            while motorX.isMoving()||motorY.isMoving()||motorZ.isMoving(),end
            toc
%             pause(.5)
            hDragonflyApp.RunSynchronizedAQ([])
        end
        
%         zSign = ~zSign;
    end
    ySign = ~ySign;
end
% xSign = ~xSign;
toc


