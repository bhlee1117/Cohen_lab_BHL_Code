function test_p = dmd_focus_test(dmd,ROI,offset,gap)
    L_lim = dmd.device.height;
    W_lim = dmd.device.width;
    L= ROI(1);
    W = ROI(2);
    
    row1 = floor(ceil(L/2)-(gap)/2)+offset(2);
    row2 = ceil(ceil(L/2)+(gap+1)/2)+offset(2);
    col1 = W/4+offset(1);
    col2 = W*3/4+offset(1);
    
test_p = false(L_lim,W_lim);
test_p([row1 row2],col1:col2)=true;
figure;imshow(test_p,[])

end