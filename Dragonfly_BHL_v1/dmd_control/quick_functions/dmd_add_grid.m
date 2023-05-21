function dmd_add_grid(app,spacing,offset)
if ~isfield(app.function_data,'grid_line')|| ~isvalid(app.function_data.grid_line)
    figure(app.fig_mask_draw)
    hold on;
    app.function_data.grid_line = line(nan,nan,'color','r');
end    

ref_im = app.ref_im;

[nrows,ncols] = size(ref_im);

grid_x = offset(1):spacing:ncols;
grid_y = offset(2):spacing:nrows;

plt_dat_x = [repmat([1; ncols],1,length(grid_y)) ...
            [grid_x;grid_x]];
plt_dat_y = [[grid_y;grid_y]  ...
            repmat([1; nrows],1,length(grid_x))];


 app.function_data.grid_line.XData = reshape([plt_dat_x;nan(1,size(plt_dat_x,2))],[],1);
 app.function_data.grid_line.YData = reshape([plt_dat_y;nan(1,size(plt_dat_y,2))],[],1);

