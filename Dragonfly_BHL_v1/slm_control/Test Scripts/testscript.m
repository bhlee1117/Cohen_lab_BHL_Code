        tic
        znew = commonz;
        znew(it_coeff) = znew(it_coeff) + zscale_init(it_coeff).*zsteps(it_steps);
        
        cplxit = zeroimg;
        for it = 1:size(ttd,2)
            phit = zeroimg;
            phit = phit + mean(zpolys124*ttd(:,it));
            phit = phit + mean(zpolys1end*znew);
            cplxit = cplxit + exp(1i*(phit));
        end

        slmph1 = angle(cplxit);
        slmph1 = slmph1*128/pi+128;
        slmph1 = max(min(slmph1,255),0);
        slmph1 = uint8(slmph1);
        toc