function result=run_CaImAn(Y,K,SNR,cell_list,argz,rangez)

Y = Y - min(Y(:)); 
if ~isa(Y,'single');    Y = single(Y);  end         % convert to single

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;             
%total number of pixels
gcp;
tau = 10;                                          % std of gaussian kernel (half size of neuron) 
pp = 2;

options = CNMFSetParms(...   
    'd1',d1,'d2',d2,...                         % dimensionality of the FOV
    'p',pp,...                                   % order of AR dynamics    
    'gSig',tau,...                              % half size of neuron
    'merge_thr',0.85,...                        % merging threshold  
    'nb',2,...                                  % number of background components    
    'min_SNR',SNR,...                             % minimum SNR threshold
    'space_thresh',1,...                      % space correlation threshold
    'cnn_thr',0.5,...            % threshold for CNN classifier    
    'deconv_method','constrained_foopsi'...
    );                                  %'deconv_method','project'...
% Data pre-processing

[P,Y] = preprocess_data(Y,pp);
% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

% display centers of found components
Cn =  correlation_image(Y); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;


    % manually refine components (optional)
refine_components = false ;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
[match_list unmatched_list]=coordinate_match(cell_list(:,[2 1 3]),center,15,argz,rangez,[]);
add_cell_list=cell_list(unmatched_list{1},[2 1 3]);
add_cell_list(find(add_cell_list(:,3)<argz-rangez | add_cell_list(:,3)>argz+rangez),:)=[];
[Ain,Cin,center,CC] = refine_components_centroid_BH(Y,Ain,Cin,center,Cn,tau,options,unmatched_list{2},add_cell_list(:,1:2)/2);
% imagesc(mean(Y,3))
% axis equal tight
% hold all
% plot(center(:,2),center(:,1),'marker','o','markersize',6,'color',[1 0 0],'linestyle','none');
% update spatial components
Yr = reshape(Y,d,T);
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

% classify components

rval_space = classify_comp_corr(Y,A,C,b,f,options);
ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test
                                        % this test will keep processes
                                        
% further classification with cnn_classifier
try  % matlab 2017b or later is needed
    [ind_cnn,value] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
end     
                            
% event exceptionality

fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
ind_exc = (fitness < options.min_fitness);

% select components

keep = (ind_corr | ind_cnn) & ind_exc;

% display kept and discarded components
A_keep = A(:,keep);
C_keep = C(keep,:);
figure;
    subplot(121); montage(extract_patch(A(:,keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15]);
        title('Kept Components');
    subplot(122); montage(extract_patch(A(:,~keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15])
        title('Discarded Components');
% merge found components
[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A_keep,b,C_keep,f,P,S,options);

%
display_merging = 1; % flag for displaying merging example
if and(display_merging, ~isempty(merged_ROIs))
    i = 1; %randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(A_keep(:,merged_ROIs{i}(j)),d1,d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(C_keep(merged_ROIs{i},:),[],2))\C_keep(merged_ROIs{i},:))'); 
            hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

% refine estimates excluding rejected components

Pm.p = 2;    % restore AR v alue
[A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
[C2,f2,P2,S2,YrA] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);


% do some plotting

[A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)

figure;
[Coor,json_file] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

% display components

%plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options);
result.C_df=C_df;
result.A_or=A_or;
result.C_or=C_or;
result.S_or=S_or;
result.P_or=P_or;
result.b2=b2;
result.f2=f2;
result.Cn=Cn;
result.Coor=Coor;
%result.Yr=Yr;

% make movie
options.save_avi=0;
 options.make_avi=0;
 result.options=options;
 %save([pth fnm(1:end-4) '_rst.mat'],'A_or','C_or','b2','f2','Yr','Coor','options','json_file','-v7.3');
 %save([save_pth fnm(1:end-4) '_rst_compressed.mat'],'C_df','A_or','C_or','b2','f2','Cn','Coor','options','json_file','-v7.3');
if (0)  
    make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)
end
