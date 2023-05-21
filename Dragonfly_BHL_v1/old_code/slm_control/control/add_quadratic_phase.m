function add_quadratic_phase(app)

x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);

lambda = .594e-3;
k=2*pi/lambda;
f=app.fmmEditField.Value; %mm

out = mod(double(app.current_mask)+ k/2/f*(X.^2+Y.^2),2*pi); % lens

app.SLM.project(out);