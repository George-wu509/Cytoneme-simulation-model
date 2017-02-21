function run_cytoneme_folder2(key,no)
% -------------------------------------------------------------------------
% updated 11/16 2015 George
% -------------------------------------------------------------------------

% --[manage folder name]---------------------------------------------------

% ==== [Parameters setting]================================================
pp.amount=[1 2 3];       % transport morphogen through cytoneme, defalut=1.5 
pp.theta=1;          % p.theta=1~12, defalut=1 
pp.length=[20 40];        % length of const cytoneme value, defalut=20 
pp.prob_format=[0.5 0.8];  % probability of cytoneme formation, defalut=0.2 
pp.gradient_L=[10];    % cytoneme length gradient coefficient, defalut=10
pp.gradient_P2=5;   % cytoneme probability gradient coefficient, defalut=20 
pp.gradient_P3=10;   % cytoneme probability gradient coefficient, defalut=10 
pp.gradient_P4=[50];   % cytoneme probability gradient coefficient, defalut=10 
pp.dfun_option=2;    % d-effect options, defalut=2 
pp.amount_dist{1}=[0,5,10,20;8,4,2,1]; % cutoneme number depend on distance, defalut=[0,5,10,20;10,5,2,1]
pp.max_gener=10;     % max generation, defalut=50 
pp.sta_save_period=5;     % figures and record every generations
pp.nn=[10];            % save figures generation in sta_all
pp.pitch=1;          % script 3 repeat times, defalut=1
pp.r=0.67;           % cell radius, defalut=0.67 
pp.area_rx=40;        % search space x length, defalut=5 
pp.area_ry=0;        % search space x length, defalut=0

% add parameters (051112)
pp.CA=1;
pp.stop_format=[0.2 0.8];  % probability of cytoneme formation, defalut=0.2
pp.gradient_P5=[1 50];   % cytoneme probability gradient coefficient, defalut=10
pp.capping=1;

%pp.Copy_list={'1_1_1_1'};
pp.Copy_list={'1_1_1_1','1_1_1_2'};

pp.name='[RUN16]';
% PA()	PL()	TH()	MA()
% 1 (const)	1 (const)	1 (const)	(1,2)
% p.amount,   p.length,   p.prob_format   p.stop_format, P5   (p.amount_dist)
% ==== [Parameters setting]================================================
[p_pre,add_folder_name]=compare_pp(pp);

pp.create_fun_file={'create3Dhex_dir.m','figdppcyto.m','cyto_position.m'};
pp.im_dig=300;
pp.im_fo1='-dtiff';
pp.im_fo2='tif';
pp.tf=ispc;
pp.nn=[10];            % save figures generation in sta_all
%% ----------------------------------------------------------------------------------------------------------

if nargin==0
    run_function=1;    % run simulation
    run_figures=0;     % draw figures
    create3Dhex_no=1;  % choice of create3Dhex_dir(N)
else
    if key==1
        run_function=1;    % run simulation
        run_figures=0;     % draw figures
        create3Dhex_no=no;  % choice of create3Dhex_dir(N)
    elseif key==2
        run_function=0;    % run simulation
        run_figures=1;     % draw figures
        create3Dhex_no=no;  % choice of create3Dhex_dir(N)
    elseif key==3
        run_function=1;    % run simulation
        run_figures=1;     % draw figures
        create3Dhex_no=no;  % choice of create3Dhex_dir(N)
    else
    end
end
pp.dpp_y=[];
pp.lambda=[];
pp.cmean=[];
pp.name='';
dirinfo = dir();pp.dirinfo=dirinfo;
dirinfo(~[dirinfo.isdir]) = [];
folder=pwd;pp.folder=folder;pp.add_folder_name=add_folder_name;
if run_function==1
    newfolder(pp,create3Dhex_no,folder);
end
if exist('Fig','dir')==0
   mkdir('Fig'); 
end
if pp.tf==1
    pp.fig_folder=[folder '\Fig'];
else
    pp.fig_folder=[folder '/Fig'];
end


dirinfo = dir();dirinfo(~[dirinfo.isdir]) = [];
for i=1:size(dirinfo)
    thisdir = dirinfo(i).name;
    if pp.tf==1
        chk_run=dir(fullfile(folder,'\',thisdir));[c,d]=size(chk_run);
    else
        chk_run=dir(fullfile(folder,'/',thisdir));[c,d]=size(chk_run);
    end
if thisdir(1,end)=='.'||strcmp(thisdir,'Fig')==1
else
    %try
    fprintf('RUN: %s\n', thisdir);
        eval(['cd(''' thisdir ''')']);
        if run_function==1&&exist('cytoneme.mat','file')==0
            if exist('n.mat')==0
                p_pre=add_pre(p_pre,thisdir);
                create3Dhex_dir(create3Dhex_no,p_pre);
            else
                load('n.mat');
                p_prejj=add_pre(p_pre{jj},thisdir);
                create3Dhex_dir(create3Dhex_no,p_prejj);
                delete('n.mat');
            end
        end
        if run_figures==1
            A=strfind(thisdir,'_');A=[A size(thisdir,2)];
            figure_title=thisdir(1:A(1)-1);
            for k=2:size(A,2)
                if pp.tf==1
                    figure_title=[figure_title '\' thisdir(A(k-1):A(k)-1)];
                else
                    figure_title=[figure_title '/' thisdir(A(k-1):A(k)-1)];
                end
            end
            figure_title=[figure_title ')'];
            if ~exist('Figures', 'dir')
            mkdir('Figures');
            end
            if exist('cytoneme.mat')~=0
            pp=createfigure(1,figure_title,pp,thisdir);%createfigure(2,figure_title,pp,thisdir);
            pp=createfigure4(1,figure_title,pp);%createfigure4(2,figure_title,pp);
            pp=createfigure6(0,1,figure_title,pp);%createfigure6(0,2,figure_title,pp);
            for nn=pp.nn
                %createfigure6(nn,1,figure_title,pp);createfigure6(nn,2,figure_title,pp);
                pp=createfigure2(nn,figure_title,pp);
                %createfigure3(nn,1,figure_title,pp);createfigure3(nn,2,figure_title,pp);
                %figdppcyto(nn,figure_title,pp);
            end
            pp=createfigure5(1,figure_title,pp);
            %show_figure(5);
            else
                display('no cytoneme.mat!');
            end
        end
        eval(['cd(''' folder ''')']);
end

end
eval(['cd(''' pp.fig_folder ''')']);
save('Summ.mat','pp');
eval(['cd(''' folder ''')']);
a=strfind(folder,'[RUN');
if isempty(a)~=1
    root_folder=folder(1:a-1);eval(['cd(''' root_folder ''')']);
    movefile(folder(a:end),[folder(a:end) 'ok'],'f');
    eval(['display(''Finish ' folder(a:end) ''');']);
end
end
function pp=createfigure(m,thisdir,pp,thisdir2)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 13-Aug-2015 15:34:05
load cytoneme.mat
YMatrix1=[];
folder_o=pwd;
n=size(sta_all,2);nnn=0;
for i=1:n
    if m==1&&mod(i,pp.sta_save_period)==0
        YMatrix1=[YMatrix1;sta_all{1,i}.dpp_y];nnn=nnn+1;
    elseif m==2&&mod(i,pp.sta_save_period)==0
        YMatrix1=[YMatrix1;sta_all{1,i}.cyto_y];nnn=nnn+1;
    end
end
n=nnn;
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.0828313253012048 0.10230179028133 0.876506024096386 0.82269820971867]);
box(axes1,'on');
hold(axes1,'all');

% Create xlabel
title(thisdir);
xlabel('Cell # (anterior-posterior axis)');

% Create ylabel
if m==1
    ylabel('[Dpp]');
elseif m==2
    ylabel('Cyto#');
end

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1','Parent',axes1,'LineStyle',':');

if n>1
    for k=1:n-1
        eval(['set(plot1(' num2str(k) '),''DisplayName'',''time= ' num2str(k*pp.sta_save_period) ''');']);
    end
    eval(['set(plot1(' num2str(n) '),''Marker'',''x'',''DisplayName'',''time=' num2str(n*pp.sta_save_period) ''',''LineStyle'',''-'');']);
else
    eval(['set(plot1(' num2str(n) '),''Marker'',''x'',''DisplayName'',''time=' num2str(n*pp.sta_save_period) ''',''LineStyle'',''-'');']);
end

% Create legend
legend1 = legend(axes1,'show');
%set(legend1,...
%    'Position',[0.758047289830463 0.643037343761981 0.176174496644295 0.21304347826087]);
cd Figures;
if m==1
    saveas(figure1,'fig_Dpp_time12345','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_Dpp_time12345' '.' pp.im_fo2 ''');']);
elseif m==2
    saveas(figure1,'fig_cyto_time12345','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_time12345' '.' pp.im_fo2 ''');']);
end
eval(['cd(''' pp.fig_folder ''')']);
saveas(figure1,thisdir2,'fig');
eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' thisdir2 '.' pp.im_fo2 ''');']);
close(figure1);
eval(['cd(''' folder_o ''')']);
end % [Dpp], cyto dynamically figures
function pp=createfigure2(i,thisdir,pp)
%CREATEFIGURE(CDATA1)
load cytoneme.mat
sta=check_cyto(sta_all{1,i},p);
cdata1=sta.cyto_matrix;
if isempty(cdata1)~=1
folder_o=pwd;

figure2 = figure;

% Create axes
axes1 = axes('Parent',figure2,'YDir','reverse',...
    'Position',[0.0421686746987952 0.0793650793650794 0.91566265060241 0.873015873015872],...
    'Layer','top');
%% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0 size(cdata1,2)]);
%% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 size(cdata1,1)]);
box(axes1,'on');
hold(axes1,'all');
title([thisdir 'Time = ' num2str(i)]);
% Create image
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
xlabel('Cell # (anterior-posterior axis)');
ylabel('cytoneme no.');
cd Figures;
file_name=['cyto_map_' num2str(i)];
saveas(figure2,file_name,'fig');
%eval(['print(figure2,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' file_name '.' pp.im_fo2 ''');']);
close(figure2);
eval(['cd(''' folder_o ''')']);
else
end
end % cytoneme distribution map
function pp=createfigure3(i,m,thisdir,pp)
%CREATEFIGURE(STA_SAVE1, STA_SAVE2, STA_SAVE3)
load cytoneme.mat
if m==1
    sta_save1=sta_all{1,i}.cyto_length_mean;
    sta_save2=sta_all{1,i}.cyto_length_std;
elseif m==2
    sta_save1=sta_all{1,i}.cyto_theta_mean;
    sta_save2=sta_all{1,i}.cyto_theta_std;
end

folder_o=pwd;
figure3 = figure;

% Create axes
axes1 = axes('Parent',figure3,...
    'Position',[0.102521008403361 0.11 0.845378151260505 0.815]);
box(axes1,'on');
hold(axes1,'all');

% Create title
title([thisdir 'Time = ' num2str(i)]);

% Create xlabel
xlabel('Cell # (anterior-posterior axis)');

% Create ylabel
ylabel('Cytoneme mean length');

% Create errorbar
if m==1
errorbar(sta_save1,sta_save2,...
    'DisplayName','sta_save{1,5}.cyto_length_std vs sta_save{1,5}.cyto_length_mean',...
    'Color',[0 0.498039215803146 0]);
plot(sta_save1,'LineWidth',2,'Color',[0 0 1],...
    'DisplayName','sta_save{1,5}.cyto_length_mean');
file_name=['cyto_length_' num2str(i)];
elseif m==2
    errorbar(sta_save1,sta_save2,...
    'DisplayName','sta_save{1,5}.cyto_length_std vs sta_save{1,5}.cyto_length_mean',...
    'Color',[0 0.498039215803146 0]);
plot(sta_save1,'LineWidth',2,'Color',[0 0 1],...
    'DisplayName','sta_save{1,5}.cyto_length_mean');
file_name=['cyto_angle_' num2str(i)];    
end

cd Figures;
saveas(figure3,file_name,'fig');
%eval(['print(figure3,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' file_name '.' pp.im_fo2 ''');']);

close(figure3);
eval(['cd(''' folder_o ''')']);
end % cytoneme angle and legnth vs x figures
function pp=createfigure4(m,thisdir,pp)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 13-Aug-2015 15:34:05
load cytoneme.mat
YMatrix1=[];
folder_o=pwd;
n=size(sta_all,2);%dpp_max=zeros(1,n+1);
nnn=0;
for i=1:n
    if m==1&&mod(i,pp.sta_save_period)==0
        A1=sta_all{1,i}.dpp_y(1,60:end-5);
        dpp_max(1,i)=max(A1);A2=A1/dpp_max(1,i);
        YMatrix1=[YMatrix1;A2];nnn=nnn+1;
    elseif m==2&&mod(i,pp.sta_save_period)==0
        A1=sta_all{1,i}.dpp_y(1,60:end-5);
        dpp_max(1,i)=max(A1);YMatrix1=[YMatrix1;A1];nnn=nnn+1;
    end
end
n=nnn;
if m==2
    dpp_max(1,n+1)=max(dpp_max(1,1:n));
    YMatrix1=YMatrix1./dpp_max(1,n+1);
end
a=size(A1,2);
XMatrix1=0:0.67*sqrt(3):0.67*sqrt(3)*(a-1);
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.0828313253012048 0.10230179028133 0.876506024096386 0.82269820971867]);
box(axes1,'on');
hold(axes1,'all');

% Create xlabel
title(thisdir);
xlabel('distance to source [um]');

% Create ylabel
if m==1
    ylabel('normalized [Dpp]');
elseif m==2
    ylabel('normalized [Dpp]');
end

% Create multiple lines using matrix input to plot
plot1 = plot(XMatrix1',YMatrix1','Parent',axes1,'LineStyle',':');

if n>1
    for k=1:n-1
        eval(['set(plot1(' num2str(k) '),''DisplayName'',''time= ' num2str(k*pp.sta_save_period) ''');']);
    end
    eval(['set(plot1(' num2str(n) '),''Marker'',''x'',''DisplayName'',''time=' num2str(n*pp.sta_save_period) ''',''LineStyle'',''-'');']);
else
    eval(['set(plot1(' num2str(n) '),''Marker'',''x'',''DisplayName'',''time=' num2str(n*pp.sta_save_period) ''',''LineStyle'',''-'');']);
end



% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.758047289830463 0.643037343761981 0.176174496644295 0.21304347826087]);
cd Figures;
if m==1
    saveas(figure1,'fig_Dpp_time12345','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_Dpp_normal' '.' pp.im_fo2 ''');']);
elseif m==2
    saveas(figure1,'fig_cyto_time12345','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_normal_all' '.' pp.im_fo2 ''');']);
end
close(figure1);
eval(['cd(''' folder_o ''')']);
end % normalized Dpp half figure 
function pp=createfigure5(m,thisdir,pp)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 13-Aug-2015 15:34:05
load cytoneme.mat
YMatrix1=[];
folder_o=pwd;
n=size(sta_all,2);
for i=n
    if m==1 %&&mod(i,pp.sta_save_period)==0
        A1=sta_all{1,i}.dpp_y(1,60:end-5);
        dpp_max=max(A1);A2=A1/dpp_max;
        YMatrix1=[YMatrix1;A2];
    elseif m==2 %&&mod(i,pp.sta_save_period)==0
        A1=sta_all{1,i}.dpp_y(1,60:end-5);
        dpp_max(1,i)=max(A1);YMatrix1=[YMatrix1;A1];
    end
end
if m==2
    dpp_max(1,n+1)=max(dpp_max(1,1:n));
    YMatrix1=YMatrix1./dpp_max(1,n+1);
end
a=size(A1,2);
XMatrix1=0:0.67*sqrt(3):0.67*sqrt(3)*(a-1);
[f,gof] = fit(XMatrix1',YMatrix1','exp1');
sta_all{1,i}.f=f;
sta_all{1,i}.lambda=-1/f.b;
sta_all{1,i}.gof=gof;
add_variable(sta_all);

text1=strfind(thisdir,'(');text2=strfind(thisdir,')');
legend1_name=thisdir(text1(1)+1:text2(end)-1);

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.0828313253012048 0.10230179028133 0.876506024096386 0.82269820971867]);
box(axes1,'on');
hold(axes1,'all');

% Create xlabel
%thisdir=[thisdir ', decay length=' num2str(-1/f.b) ' (um)'];
title(thisdir);
%xlabel('Cell # (anterior-posterior axis)');

% Create ylabel
%if m==1
%    ylabel('normalized intensity(each generation)');
%elseif m==2
%    ylabel('normalized intensity');
%end

% Create multiple lines using matrix input to plot
plot1 = plot(f,XMatrix1',YMatrix1');
%plot1 = plot(f,XMatrix1',YMatrix1','Parent',axes1,'LineStyle',':');
p.legend=[legend1_name ',   \lambda=' num2str(-1/f.b)];
set(plot1(1),'DisplayName',p.legend);
set(plot1(2),'DisplayName','fitted curve','LineStyle','-');
xlabel('distance to source [um]');

% Create ylabel
if m==1
    ylabel('normalized [Dpp]');
elseif m==2
    ylabel('normalized [Dpp]');
end

% Create legend
legend1 = legend(axes1,'show');
%set(legend1,...
%    'Position',[0.758047289830463 0.643037343761981 0.176174496644295 0.21304347826087]);
cd Figures;
if m==1
    saveas(figure1,'fig_Dpp_fit','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_Dpp_fit' '.' pp.im_fo2 ''');']);
elseif m==2
    saveas(figure1,'fig_cyto_time12345','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_fit' '.' pp.im_fo2 ''');']);
end
close(figure1);
eval(['cd(''' folder_o ''')']);

pp.dpp_y=[pp.dpp_y;sta_all{end}.dpp_y];
pp.lambda=[pp.lambda;sta_all{end}.lambda];
pp.f{size(pp.dpp_y,1)}=sta_all{end}.f;
pp.gof{size(pp.dpp_y,1)}=sta_all{end}.gof;
pp.name=[pp.name;{thisdir}];
end % exp fitting curve figure
function pp=createfigure6(ii,m,thisdir,pp)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 13-Aug-2015 15:34:05
load cytoneme.mat
YMatrix1=[];YMatrix2=[];
folder_o=pwd;
n=size(sta_all,2);B=zeros(n,13);nnn=0;
if ii==0
    A=1:n;
    for i=A
    if m==1&&mod(i,pp.sta_save_period)==0
        nn=size(sta_all{1,i}.cyto,2);cyto_length=zeros(1,nn);
        for k=1:nn
            cyto_length(1,k)=distance(sta_all{1,i}.cyto(3:4,k)',sta_all{1,i}.cyto(7:8,k)');
            temp=min(fix(cyto_length(1,k)/5)+1,13);
            B(i,temp)=B(i,temp)+1;
        end
        YMatrix1=[YMatrix1;B(i,:)];nnn=nnn+1;
        cyto_leng{i}=cyto_length;
    elseif m==2&&mod(i,pp.sta_save_period)==0
        nn=size(sta_all{1,i}.cyto,2);cyto_length=zeros(1,nn);
        for k=1:nn
            cyto_length(1,k)=distance(sta_all{1,i}.cyto(3:4,k)',sta_all{1,i}.cyto(7:8,k)');
            temp=min(fix(cyto_length(1,k)/5)+1,13);
            B(i,temp)=B(i,temp)+1;
        end
        YMatrix1=[YMatrix1;B(i,:)/max(B(i,:))];nnn=nnn+1;
        cyto_leng{i}=cyto_length;
    end
    end
    n=nnn;   
else
    A=ii;
    for i=A
        if m==1
            nn=size(sta_all{1,i}.cyto,2);cyto_length=zeros(1,nn);
            for k=1:nn
                cyto_length(1,k)=distance(sta_all{1,i}.cyto(3:4,k)',sta_all{1,i}.cyto(7:8,k)');
                temp=min(fix(cyto_length(1,k)/5)+1,13);
                B(i,temp)=B(i,temp)+1;
            end
            YMatrix1=[YMatrix1;B(i,:)];nnn=nnn+1;
            cyto_leng{i}=cyto_length;
        elseif m==2
            nn=size(sta_all{1,i}.cyto,2);cyto_length=zeros(1,nn);
            for k=1:nn
                cyto_length(1,k)=distance(sta_all{1,i}.cyto(3:4,k)',sta_all{1,i}.cyto(7:8,k)');
                temp=min(fix(cyto_length(1,k)/5)+1,13);
                B(i,temp)=B(i,temp)+1;
            end
            YMatrix1=[YMatrix1;B(i,:)/max(B(i,:))];nnn=nnn+1;
            cyto_leng{i}=cyto_length;
        end
    end
    n=nnn;
end
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'0-5','10-15','20-25','30-35','40-45','50-55','60-65','>65'});
box(axes1,'on');
hold(axes1,'all');

% Create xlabel
if ii==0
    title(thisdir);
    for i=A
        cmean=mean(cyto_leng{i});
        sta_all{1,i}.cmean=cmean;
    end
    add_variable(sta_all);
else
    cmean=mean(cyto_length);
    sta_all{1,ii}.cmean=cmean;
    pp.cmean=[pp.cmean;cmean];
    add_variable(sta_all);
    title([thisdir ', mean cytoneme length=' num2str(cmean) '(um)']);
end
pp.cmean=[pp.cmean;cmean];
xlabel('cytoneme length (um)');

% Create ylabel
if m==1
    ylabel('cytoneme number');
elseif m==2
    ylabel('normalized cytoneme number');
end

% Create multiple lines using matrix input to plot
if ii==0
plot1 = bar(YMatrix1','Parent',axes1,'LineStyle',':');
%set(plot1(1),'DisplayName','time=10');
%set(plot1(2),'DisplayName','time=20');
%set(plot1(3),'DisplayName','time=30');
%set(plot1(4),'DisplayName','time=40');
%set(plot1(5),'DisplayName','time=50');

for k=1:n
    eval(['set(plot1(' num2str(k) '),''DisplayName'',''time= ' num2str(k*pp.sta_save_period) ''');']);
end

% Create legend
legend1 = legend(axes1,'show');
%set(legend1,...
%    'Position',[0.758047289830463 0.643037343761981 0.176174496644295 0.21304347826087]);
cd Figures;
if m==1
    saveas(figure1,'fig_cyto_dist','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_dist' '.' pp.im_fo2 ''');']);
elseif m==2
    saveas(figure1,'fig_cyto_dist_norm','fig');
    eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_dist_norm' '.' pp.im_fo2 ''');']);
end
close(figure1);
else
    %bar(YMatrix1,'DisplayName','A');
    bar1=bar(YMatrix1);
    eval(['set(bar1,''DisplayName'',''time=' num2str(ii) ''');'])

% Create legend
legend1 = legend(axes1,'show');
%set(legend1,...
%    'Position',[0.758047289830463 0.643037343761981 0.176174496644295 0.21304347826087]);
cd Figures;
if m==1
    if ii==0
        saveas(figure1,'fig_cyto_dist','fig');
        eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_dist' '.' pp.im_fo2 ''');']);
    else
        eval(['saveas(figure1,''fig_cyto_dist_' num2str(ii) ''',''fig'');'])
        eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_dist' '_' num2str(ii) '.' pp.im_fo2 ''');']);
    end
elseif m==2
    if ii==0
        saveas(figure1,'fig_cyto_dist_norm','fig');
        eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_dist_norm' '.' pp.im_fo2 ''');']);
    else
        eval(['saveas(figure1,''fig_cyto_dist_norm_' num2str(ii) ''',''fig'');'])
        eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' 'fig_cyto_dist_norm' '_' num2str(ii) '.' pp.im_fo2 ''');']);
    end
end
close(figure1);
end
eval(['cd(''' folder_o ''')']);
end % cytoneme length distribution figure
function figdppcyto(n,thisdir,pp)

xrange=[];
yrange=[];

%  Auto-generated by MATLAB on 22-Oct-2015 04:31:57
load('cytoneme.mat');
cdata1=sta_all{1,n}.dpp;gene=n;
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YDir','reverse',...
    'Position',[0.0445168295331162 0.11 0.92399565689468 0.815],...
    'Layer','top');
%% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0.5 120.5]);
%% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0.5 50.5]);
box(axes1,'on');
hold(axes1,'all');

% Create image
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
title([thisdir 'Time = ' num2str(n)]);
% Create plot
for i=1:size(sta_all{1,gene}.cyto,2)
    X1=[sta_all{1,gene}.cyto(2,i) sta_all{1,gene}.cyto(6,i)];
    Y1=[sta_all{1,gene}.cyto(1,i) sta_all{1,gene}.cyto(5,i)];
    if isempty(xrange)~=1
        if min(xrange)<min(X1)&&max(xrange)<min(X1)
        elseif min(xrange)>max(X1)&&max(xrange)>max(X1)
        else
            if isempty(yrange)~=1
                if min(yrange)<min(Y1)&&max(yrange)<min(Y1)
                elseif min(yrange)>max(Y1)&&max(yrange)>max(Y1)
                else
                    eval(['plot(X1,Y1,''Visible'',''on'',''DisplayName'',''(' num2str(Y1(1)) ',' num2str(X1(1)) ') to (' num2str(Y1(2)) ',' num2str(X1(2)) ')'');']);
                end                
            else
                eval(['plot(X1,Y1,''Visible'',''on'',''DisplayName'',''(' num2str(Y1(1)) ',' num2str(X1(1)) ') to (' num2str(Y1(2)) ',' num2str(X1(2)) ')'');']);
            end
        end
    else
        if isempty(yrange)~=1
            if min(yrange)<min(Y1)&&max(yrange)<min(Y1)
            elseif min(yrange)>max(Y1)&&max(yrange)>max(Y1)
            else
                eval(['plot(X1,Y1,''Visible'',''on'',''DisplayName'',''(' num2str(Y1(1)) ',' num2str(X1(1)) ') to (' num2str(Y1(2)) ',' num2str(X1(2)) ')'');']);
            end                
        else
            eval(['plot(X1,Y1,''Visible'',''on'',''DisplayName'',''(' num2str(Y1(1)) ',' num2str(X1(1)) ') to (' num2str(Y1(2)) ',' num2str(X1(2)) ')'');']);
        end
    end
end
% save figures
cd Figures;
file_name=['cyto_location_' num2str(n)];
saveas(figure1,file_name,'fig');
eval(['print(figure1,''-r' num2str(pp.im_dig) ''',''' pp.im_fo1 ''',''' file_name '.' pp.im_fo2 ''');']);
close(figure1);
eval(['cd(''' folder_o ''')']);
end
function newfolder(pp,create3Dhex_no,folder)
n=size(pp.Copy_list,2)*size(pp.add_folder_name,2);
for i=1:size(pp.Copy_list,2)
    for jj=1:size(pp.add_folder_name,2)
        if create3Dhex_no==4
            aa='(T)';
        elseif create3Dhex_no==6
            aa='(M)';
        elseif create3Dhex_no==7
            aa='(TM)';
        else
            aa='';
        end
        newfolder=[pp.Copy_list{1,i} aa pp.add_folder_name{1,jj}];
        chk1=0;
        for j=1:size(pp.dirinfo,1)
            if  strcmp(newfolder,pp.dirinfo(j,1).name)==1
                chk1=1; %existed folder
            end
        end
        if chk1==1
            if newfolder(end)==')'&&newfolder(end-2)=='('
                N=str2num(A(end-1));
                newfolder=[newfolder '(' num2str(N+1) ')'];
            else
                newfolder=[newfolder '(' num2str(1) ')'];
            end
        %else
        end
        eval(['mkdir(''' newfolder ''');']);
        if pp.tf==1
            for jjj=1:size(pp.create_fun_file,2)
                copyfile(pp.create_fun_file{jjj},newfolder);dd_folder=[folder '\' newfolder '\' 'n.mat'];
            end
        else     
            for jjj=1:size(pp.create_fun_file,2)
                copyfile(pp.create_fun_file{jjj},newfolder);dd_folder=[folder '/' newfolder '/' 'n.mat'];
            end
        end
            save(dd_folder,'jj');
            eval(['display(''create ' newfolder ''');']);
        %end
        eval(['cd(''' pp.folder ''');']);
    end
end
end
function edit_file(create_fun_file,Copy_list)
create_fun_file='testtest.m';
if exist(create_fun_file)~=0
fid=fopen(create_fun_file,'r');
y=0;
while feof(fid)==0
    tline=fgetl(fid);
    matches=findstr(tline,'p.AN=');
    num=length(matches);
    if num>0
        fprintf(fid, '%s', 'test');
    end
end
end
end
function p_pre=add_pre(p_pre,thisdir)
a=size(p_pre,2);
eval(['p_pre{a+1}=''p.AN=' thisdir(1) ';'';'])
eval(['p_pre{a+2}=''p.LE=' thisdir(3) ';'';'])
eval(['p_pre{a+3}=''p.PF=' thisdir(5) ';'';'])
eval(['p_pre{a+4}=''p.MA=' thisdir(7) ';'';'])
end
function [p_pre,add_folder_name]=compare_pp(pp)
dp.amount=1.5;       % transport morphogen through cytoneme, defalut=1.5 
dp.theta=1;          % p.theta=1~12, defalut=1 
dp.length=40;        % length of const cytoneme value, defalut=20 
dp.prob_format=0.2;  % probability of cytoneme formation, defalut=0.2 
dp.gradient_L=10;    % cytoneme length gradient coefficient, defalut=10
dp.gradient_P2=20;   % cytoneme probability gradient coefficient, defalut=20 
dp.gradient_P3=10;   % cytoneme probability gradient coefficient, defalut=10 
dp.gradient_P4=50;   % cytoneme probability gradient coefficient, defalut=10 
dp.dfun_option=2;    % d-effect options, defalut=2 
dp.amount_dist{1}=[0,5,10,20;10,5,2,1]; % cutoneme number depend on distance, defalut=[0,5,10,20;10,5,2,1]
dp.max_gener=50;     % max generation, defalut=50
pp.sta_save_period=10;     % figures and record every generations
dp.pitch=1;          % script 3 repeat times, defalut=1
dp.r=0.67;           % cell radius, defalut=0.67 
dp.area_rx=40;        % search space x length, defalut=5 
dp.area_ry=0;        % search space x length, defalut=0
dp.CA=1;
dp.stop_format=[0.2];  % probability of cytoneme formation, defalut=0.2
dp.gradient_P5=[50];   % cytoneme probability gradient coefficient, defalut=10
dp.capping=0;

sizeDp = length(fieldnames(dp));sizePp = length(fieldnames(pp));
k1=1;add_folder_name=[];p_name=fieldnames(dp);p_term=[];
if sizePp>=sizeDp    
for i=1:sizeDp
    eval(['a1=dp.' p_name{i} ';']);eval(['a2=pp.' p_name{i} ';']);
    if strcmp(p_name{i},'amount_dist')==0
    if size(a2,2)>1||isempty(find(a1==a2, 1))
        p_term{k1,1}=p_name{i};p_term{k1,2}=a2;p_term{k1,3}=size(a2,2);
        k1=k1+1;
    end
    else
        if iscell(a2)==1
            if size(a2,2)==1
                if isempty(find(a2{1}-a1{1}, 1))==0
                    p_term2=a2;p_term2_s=size(p_term2,2);
                else
                    p_term2=0;p_term2_s=0;
                end
            else
                p_term2=a2;p_term2_s=size(p_term2,2);
            end
        elseif isempty(find(a2-a1{1}, 1))==0
            p_term2{1}=a2;p_term2_s=1;
        else
            p_term2=0;p_term2_s=0;
        end
    end
end

% calculate s, s_list
if isempty(p_term)==1 % no other parameters control
    if p_term2_s==0
        s=p_term2_s;s_list=0;
    else
        s=p_term2_s;s_list=1:p_term2_s;
    end
else
    if p_term2_s==0  % other parameters control, no amount_dist.
        s=1;
        for ii=1:size(p_term,1)
            s=s.*p_term{ii,3};
        end
        s_list=zeros(size(p_term,1),s);
        s_list(1,:)=rem(0:s-1,p_term{1,3})+1;
        if size(p_term,1)>1
            cc=0;
        for ii=2:size(p_term,1)
            num=1;s_list(ii,1)=num;
            if p_term{ii,3}==1
                s_list(ii,:)=1;cc=cc+1;
            else           
            for kk=2:s
                if s_list(ii-1-cc,kk-1)==s_list(ii-1-cc,end)&&s_list(ii-1-cc,kk)==1
                    num=num+1;
                    if num>p_term{ii,3}
                        num=1;
                    end
                end
                s_list(ii,kk)=num;
            end
            cc=0;
            end
        end
        end
    else   % other parameters control, with amount_dist.
        s=p_term2_s;
        for ii=1:size(p_term,1)
        s=s.*p_term{ii,3};
        end
            s_list=zeros(size(p_term,1)+1,s);
            s_list(1,:)=rem(0:s-1,p_term2_s)+1;
            
            if size(p_term,1)>=1
                cc=0;
            for ii=2:size(p_term,1)+1
                num=1;s_list(ii,1)=num;
                if p_term{ii-1,3}==1
                    s_list(ii,:)=1;cc=cc+1;
                else
                for kk=2:s
                    if s_list(ii-1-cc,kk-1)==s_list(ii-1-cc,end)&&s_list(ii-1-cc,kk)==1
                        num=num+1;
                        if num>p_term{ii-1,3}
                            num=1;
                        end
                    end
                    s_list(ii,kk)=num;
                end
                cc=0;
                end
            end
            end
    end    
end

% output p_pre, add_folder_name
if s~=0
for i=1:s
    add_folder_name{i}='';
    if p_term2_s==0
        for j=1:size(s_list,1)
            eval(['p_pre{i}{j}=''p.' p_term{j,1} '=' num2str(p_term{j,2}(s_list(j,i))) ';'';']);
            add_folder_name{i}=ad_name_fun(add_folder_name{i},p_term{j,1},p_term{j,2}(s_list(j,i)));
        end
    else
        %eval(['p_pre{i}{1}=''p.amount_dist=' num2str(s_list(1,i)) ';'';']);
        [add_folder_name{i} p_pre{i}{1}]=ad_name_fun2(add_folder_name{i},s_list(1,i),p_term2);
        for j=2:size(s_list,1)
            eval(['p_pre{i}{j}=''p.' p_term{j-1,1} '=' num2str(p_term{j-1,2}(s_list(j,i))) ';'';']);
            add_folder_name{i}=ad_name_fun(add_folder_name{i},p_term{j-1,1},p_term{j-1,2}(s_list(j,i)));
        end       
    end
end
end
end
if exist('p_pre','var')==0
    p_pre{1}='';
end
end
function add_folder_name=ad_name_fun(add_folder_name,id,num)
dpname={'m','theta','leng','PF','P1','P2','P3','dfun','a_d','maxg','repet','r','rx','ry'};
%id=p_term{j,1};num=p_term{j,2}(s_list1);
if strcmp(id,'amount')==1    
    eval(['add_folder_name=[add_folder_name ''(m=' num2str(num) ')''];']);
elseif strcmp(id,'theta')==1   
    eval(['add_folder_name=[add_folder_name ''(theta=' num2str(num) ')''];']);
elseif strcmp(id,'length')==1   
    eval(['add_folder_name=[add_folder_name ''(leng=' num2str(num) ')''];']);
elseif strcmp(id,'prob_format')==1   
    eval(['add_folder_name=[add_folder_name ''(IP=' num2str(num) ')''];']);
elseif strcmp(id,'gradient_L')==1   
    eval(['add_folder_name=[add_folder_name ''(P1=' num2str(num) ')''];']);
elseif strcmp(id,'gradient_P2')==1   
    eval(['add_folder_name=[add_folder_name ''(P2=' num2str(num) ')''];']);
elseif strcmp(id,'gradient_P3')==1   
    eval(['add_folder_name=[add_folder_name ''(P3=' num2str(num) ')''];']);
elseif strcmp(id,'gradient_P4')==1   
    eval(['add_folder_name=[add_folder_name ''(P4=' num2str(num) ')''];']);
elseif strcmp(id,'dfun_option')==1   
    eval(['add_folder_name=[add_folder_name ''(dfun=' num2str(num) ')''];']);
elseif strcmp(id,'amount_dist')==1   
    eval(['add_folder_name=[add_folder_name ''(a_d=' num2str(num) ')''];']);
elseif strcmp(id,'max_gener')==1   
    eval(['add_folder_name=[add_folder_name ''(maxg=' num2str(num) ')''];']);
elseif strcmp(id,'pitch')==1   
    eval(['add_folder_name=[add_folder_name ''(repet=' num2str(num) ')''];']);
elseif strcmp(id,'r')==1   
    eval(['add_folder_name=[add_folder_name ''(' num2str(num)/0.67 'r'')''];']);
elseif strcmp(id,'area_rx')==1   
    eval(['add_folder_name=[add_folder_name ''(rx=' num2str(num) ')''];']);
elseif strcmp(id,'area_ry')==1   
    eval(['add_folder_name=[add_folder_name ''(ry=' num2str(num) ')''];']);  
elseif strcmp(id,'CA')==1   
    eval(['add_folder_name=[add_folder_name ''(CA=' num2str(num) ')''];']);  
elseif strcmp(id,'stop_format')==1   
    eval(['add_folder_name=[add_folder_name ''(CP=' num2str(num) ')''];']); 
elseif strcmp(id,'gradient_P5')==1   
    eval(['add_folder_name=[add_folder_name ''(P5=' num2str(num) ')''];']);
elseif strcmp(id,'capping')==1   
    eval(['add_folder_name=[add_folder_name ''(cap=' num2str(num) ')''];']);
end
end
function [add_folder_name p_pre]=ad_name_fun2(add_folder_name,s_list1,p_term2)
ad=p_term2{1,s_list1};
eval(['add_folder_name=[add_folder_name ''(ad=' num2str(ad(2,4)) num2str(ad(2,3)) num2str(ad(2,2)) num2str(ad(2,1)) ')''];']);

eval(['p_pre=''p.amount_dist=[0,5,10,20;' num2str(ad(2,1)) ',' num2str(ad(2,2)) ',' num2str(ad(2,3)) ',' num2str(ad(2,4)) '] ;'';']);
end

function show_figure(nn)
load cytoneme.mat
folder_o=pwd;
p.save_every_img=1;p.show_every_img=1;
sta=sta_save{1,nn};ge=10*nn;
    if p.show_every_img==1;
        % draw mesh color distribution
        F1=figure('visible','on');[sta,p]=draw_mesh(sta,p,ge);
        % draw 3D mesh color distribution
        F2=figure('visible','on');drawsurf(sta,p,ge);
    else
        % draw mesh color distribution
        F1=figure('visible','off');[sta,p]=draw_mesh(sta,p,ge);
        % draw 3D mesh color distribution
        F2=figure('visible','off');drawsurf(sta,p,ge);
    end
    if p.save_every_img==1;
        cd Figures;
        eval(['saveas(F1,''f1_' num2str(ge) ''',''jpg'');']);
        eval(['saveas(F2,''f2_' num2str(ge) ''',''jpg'');']);
        eval(['saveas(F1,''f1_' num2str(ge) ''',''fig'');']);
        eval(['saveas(F2,''f2_' num2str(ge) ''',''fig'');']);
        eval(['cd(''' folder_o ''')']);
    end
    p=draw_cytoneme(sta,p,ge);
    %{
    if p.show_every_img==1;
        F4=figure('visible','on');
    else
        F4=figure('visible','off');
    end   
    if exist('sta_save','var')~=1
        draw_section(sta,ge);
    else
        draw_section(sta,ge,sta_save);
    end
    close(F4);  
    %}
    close(F1);close(F2);
end
function [sta,p]=draw_mesh(sta,p,ge)
% hex line
v=30:60:390;
cv=p.r*cosd(v);
sv=p.r*sind(v);
cmap=colormap(jet(p.col_range(1,2)));
dpp_maxmin(1,1)=max(max(sta.dpp));dpp_maxmin(1,2)=min(min(sta.dpp));

% 2D or 3D mesh
if p.height==0
for j=1:p.nx  % create 2D first layer
for i=1:p.ny
    if p.xc(i,j)>=0&&p.yc(i,j)>=0
    % Draw boundary line
    bl_seq=[p.xc(i,j)+cv;p.yc(i,j)+sv];
    %line(bl_seq(1,:),bl_seq(2,:),'tag','h');
    % Draw color on upper and lower surface   
    %{
    sta.rgb(i,j)=sta.dpp(i,j)+1;
    if sta.rgb(i,j)>p.col_range(1,2)
        sta.rgb(i,j)=p.col_range(1,2);
    end
    %}
    if dpp_maxmin(1,1)-dpp_maxmin(1,2)==0
        sta.rgb(i,j)=p.col_range(1,1);
    else
        sta.rgb(i,j)=round(((sta.dpp(i,j)-dpp_maxmin(1,2))/(dpp_maxmin(1,1)-dpp_maxmin(1,2)))*(p.col_range(1,2)-p.col_range(1,1))+p.col_range(1,1));
    end   
    if p.h2d_on==0;
        p.h2d(1,i,j)=patch(bl_seq(1,:)',bl_seq(2,:)',cmap(sta.rgb(i,j),1:3)); % lower surface
    else
        set(p.h2d(1,i,j),'FaceColor',cmap(sta.rgb(i,j),1:3));
    end
    end
end
end
p.h2d_on=1;
xlabel('x(A->P)');ylabel('y(V->D)');title(['dpp concentration' ' (generation ' num2str(ge) ')'])

% 3D mesh
else
for j=1:p.nx
for i=1:p.ny
    if p.xc(i,j)>=0&&p.yc(i,j)>=0
    % Draw boundary line
    bl_seq=[p.xc(i,j)+cv;p.yc(i,j)+sv];
    %line(bl_seq(1,:),bl_seq(2,:),zeros(1,size(bl_seq,2)),'tag','h');
    %line(bl_seq(1,:),bl_seq(2,:),ones(1,size(bl_seq,2))*p.height,'tag','h');
    for hi=1:size(bl_seq,2)
    %    line(ones(1,2)*bl_seq(1,hi),ones(1,2)*bl_seq(2,hi),[0,p.height],'tag','h');
    end
    % Draw color on upper and lower surface
    sta.rgb(i,j)=round(sta.dpp(i,j)+1);
    if sta.rgb(i,j)>p.col_range(1,2)
        sta.rgb(i,j)=p.col_range(1,2);
    elseif sta.rgb(i,j)<1
        sta.rgb(i,j)=1;
    end
    if p.h3d_on==0;
        p.h3d(1,(j-1)*p.ny+i)=patch(bl_seq(1,:)',bl_seq(2,:)',zeros(1,size(bl_seq,2))',cmap(sta.rgb(i,j),1:3)); % lower surface
        p.h3d(2,(j-1)*p.ny+i)=patch(bl_seq(1,:)',bl_seq(2,:)',ones(1,size(bl_seq,2))*p.height',cmap(sta.rgb(i,j),1:3)); % upper surface
        % Draw color on lateral surface 
        for k=1:size(bl_seq,2)-1
            xx=[bl_seq(1,k:k+1) bl_seq(1,k+1:-1:k)]';
            yy=[bl_seq(2,k:k+1) bl_seq(2,k+1:-1:k)]';
            p.h3d(k+2,(j-1)*p.ny+i)=patch(xx,yy,[zeros(2,1);ones(2,1)*p.height],cmap(sta.rgb(i,j),1:3)); % lateral surface
        end
    else
        set(p.h3d(1,(j-1)*p.ny+i),'FaceColor',cmap(sta.rgb(i,j),1:3));
        set(p.h3d(2,(j-1)*p.ny+i),'FaceColor',cmap(sta.rgb(i,j),1:3));
        % Draw color on lateral surface 
        for k=1:size(bl_seq,2)-1
            set(p.h3d(k+2,(j-1)*p.ny+i),'FaceColor',cmap(sta.rgb(i,j),1:3)); % lateral surface
        end
    end
    end
end
end
p.h3d_on=1;
xlabel('x(A->P) um');ylabel('y(V->D) um');zlabel('z(basal->apical) um');title(['dpp concentration' ' (generation ' num2str(ge) ')'])
end

% Display region
%colorbar('location','North');
alpha(p.alpha1);
axis equal;
%camlight right;
%axis([0 x 0 y ]);
ScreenSize=get(0,'ScreenSize');
set(gcf,'Position',[ScreenSize(3)*0.02 ScreenSize(4)*0.48 ScreenSize(3)*0.4 ScreenSize(4)*0.4]);
end
function drawsurf(sta,p,ge)
    surf(p.xindex,p.yindex,sta.dpp*0.1,sta.rgb); shading interp
    
    % Display region
    xlabel('xi(A->P)');ylabel('yi(V->D)');zlabel('[Dpp]*10^-1');title(['[Dpp] on columnar layer' ' (generation ' num2str(ge) ')' ]);
    colorbar;
    alpha(p.alpha2);
    axis equal;
    %camlight right;
    %axis([0 x 0 y ]);
    %view([-25+(i/p.sta_save_period),30]);
    view([-25,30]);
    ScreenSize=get(0,'ScreenSize');
    set(gcf,'Position',[ScreenSize(3)*0.35 ScreenSize(4)*0.1 ScreenSize(3)*0.63 ScreenSize(4)*0.6]);
end
function p=draw_cytoneme(sta,p,ge)
folder_o=pwd;
if p.show_every_img==1;
    F3=figure('visible','on');
else
    F3=figure('visible','off');
end
F3_line=findobj(F3, 'Type', 'line');
	if ~isempty(F3_line)
		delete(F3_line);
	end	
% hex line
v=30:60:390;
cv=p.r*cosd(v);
sv=p.r*sind(v);
cmap=colormap(jet(p.col_range(1,2)));

% 2D or 3D mesh
if p.height==0
for j=1:p.nx  % create 2D first layer
for i=1:p.ny
    if p.xc(i,j)>=0&&p.yc(i,j)>=0
    % Draw boundary line
    bl_seq=[p.xc(i,j)+cv;p.yc(i,j)+sv];
    %line(bl_seq(1,:),bl_seq(2,:),'tag','h');
    % Draw color on upper and lower surface
    %{
    sta.rgb(i,j)=sta.dpp(i,j)+1;
    if sta.rgb(i,j)>p.col_range(1,2)
        sta.rgb(i,j)=p.col_range(1,2);
    end
    %}
    if p.h2dh_on==0;
        p.h2dh(1,(j-1)*p.ny+i)=patch(bl_seq(1,:)',bl_seq(2,:)',cmap(sta.rgb(i,j),1:3),'EdgeColor','w'); % lower surface
    else
        set(p.h2dh(1,(j-1)*p.ny+i),'FaceColor',cmap(sta.rgb(i,j),1:3));
    end
    end
end
end
p.h2dh_on=1;

% 3D mesh
else
for j=1:p.nx
for i=1:p.ny
    if p.xc(i,j)>=0&&p.yc(i,j)>=0
    % Draw boundary line
    bl_seq=[p.xc(i,j)+cv;p.yc(i,j)+sv];
    %line(bl_seq(1,:),bl_seq(2,:),zeros(1,size(bl_seq,2)),'tag','h');
    %line(bl_seq(1,:),bl_seq(2,:),ones(1,size(bl_seq,2))*p.height,'tag','h');
    for hi=1:size(bl_seq,2)
    %    line(ones(1,2)*bl_seq(1,hi),ones(1,2)*bl_seq(2,hi),[0,p.height],'tag','h');
    end
    % Draw color on upper and lower surface
    sta.rgb(i,j)=round(sta.dpp(i,j)+1);
    if sta.rgb(i,j)>p.col_range(1,2)
        sta.rgb(i,j)=p.col_range(1,2);
    elseif sta.rgb(i,j)<1
        sta.rgb(i,j)=1;
    end
    if p.h3dh_on==0;
        p.h3dh(1,(j-1)*p.ny+i)=patch(bl_seq(1,:)',bl_seq(2,:)',zeros(1,size(bl_seq,2))',cmap(sta.rgb(i,j),1:3),'EdgeColor','w'); % lower surface
        p.h3dh(2,(j-1)*p.ny+i)=patch(bl_seq(1,:)',bl_seq(2,:)',ones(1,size(bl_seq,2))*p.height',cmap(sta.rgb(i,j),1:3),'EdgeColor','w'); % upper surface
        % Draw color on lateral surface 
        for k=1:size(bl_seq,2)-1
            xx=[bl_seq(1,k:k+1) bl_seq(1,k+1:-1:k)]';
            yy=[bl_seq(2,k:k+1) bl_seq(2,k+1:-1:k)]';
            p.h3dh(k+2,(j-1)*p.ny+i)=patch(xx,yy,[zeros(2,1);ones(2,1)*p.height],cmap(sta.rgb(i,j),1:3),'EdgeColor','w'); % lateral surface
        end
    else
        set(p.h3dh(1,(j-1)*p.ny+i),'FaceColor',cmap(sta.rgb(i,j),1:3));
        set(p.h3dh(2,(j-1)*p.ny+i),'FaceColor',cmap(sta.rgb(i,j),1:3));
        % Draw color on lateral surface 
        for k=1:size(bl_seq,2)-1
            set(p.h3dh(k+2,(j-1)*p.ny+i),'FaceColor',cmap(sta.rgb(i,j),1:3)); % lateral surface
        end
    end
    end
end
end
p.h3dh_on=1;
end

% draw cytoneme
if isempty(sta.cyto)
else
movelist=[0;0;0];
cy_n=size(sta.cyto,2);
for c=1:cy_n
    [d1,d2]=min(sum(abs(movelist(1:2,:)-sta.cyto(1:2,c)*ones(1,size(movelist,2)))));
    [d3,d4]=min(sum(abs(movelist(1:2,:)-sta.cyto(5:6,c)*ones(1,size(movelist,2)))));
    if d1==0
        movelist(3,d2)=movelist(3,d2)+1;
    else
        movelist=[movelist [sta.cyto(1:2,c);0]];
    end
    if d3==0
        movelist(3,d4)=movelist(3,d4)+1;
    else
        movelist=[movelist [sta.cyto(5:6,c);0]];
    end
    dd=0.2;
    if p.height==0
        line([sta.cyto(3,c)+movelist(3,d2)*dd sta.cyto(7,c)+movelist(3,d4)*dd],[sta.cyto(4,c)+movelist(3,d2)*dd sta.cyto(8,c)+movelist(3,d4)*dd],'LineWidth',1.0);
    else
        line([sta.cyto(3,c)+movelist(3,d2)*dd sta.cyto(7,c)+movelist(3,d4)*dd],[sta.cyto(4,c)+movelist(3,d2)*dd sta.cyto(8,c)+movelist(3,d4)*dd],[p.height p.height],'LineWidth',1.0);
    end
end
end

% Display region
%colorbar;
xlabel('x(A->P) um');ylabel('y(V->D) um');zlabel('z(basal->apical) um');title(['concentration events' ' (generation ' num2str(ge) ')'])
view([-20,50]);
alpha(0.1);
axis equal;
%camlight right;
%axis([0 x 0 y ]);
ScreenSize=get(0,'ScreenSize');
set(gcf,'Position',[ScreenSize(3)*0.55 ScreenSize(4)*0.5 ScreenSize(3)*0.42 ScreenSize(4)*0.38]);
    if p.save_every_img==1;
        cd Figures;
        eval(['saveas(F3,''f3_' num2str(ge) ''',''jpg'');']);
        eval(['saveas(F3,''f3_' num2str(ge) ''',''fig'');']);
        eval(['cd(''' folder_o ''')']);
    end
    close(F3);
end
function draw_section(sta,ge,sta_save)

if exist('sta_save','var')==1
    if isempty(sta_save)==1
        A1(1,:)=mean(sta.dpp,1);
        A2(:,1)=mean(sta.dpp,2);
    else
        for j=1:size(sta_save,2)
        A1(j,:)=mean(sta_save{1,j}.dpp,1);
        A2(:,j)=mean(sta_save{1,j}.dpp,2);
        end
    end
else
    A1(1,:)=mean(sta.dpp,1);
    A2(:,1)=mean(sta.dpp,2);
end
if size(A1,1)>1
    subplot(2,1,1);plot(A1(1:end-1,:)','LineStyle',':');hold on;plot(A1(end,:)','LineStyle','-','LineWidth',1);hold off;
    title(['mean Dpp concentration through AP axis' ' (generation ' num2str(ge) ')']);
    xlabel('x(P->A)');ylabel('mean [Dpp]');
    subplot(2,1,2);plot(A2(:,1:end-1),'LineStyle',':');hold on;plot(A2(:,end),'LineStyle','-','LineWidth',1);hold off;
    title(['mean Dpp concentration through DV axis' ' (generation ' num2str(ge) ')']);
    xlabel('y(V->D)');ylabel('mean [Dpp]');
else
    subplot(2,1,1);plot(A1');
    title(['mean Dpp concentration through AP axis' ' (generation ' num2str(ge) ')']);
    xlabel('x(P->A)');ylabel('mean [Dpp]');
    subplot(2,1,2);plot(A2);
    title(['mean Dpp concentration through DV axis' ' (generation ' num2str(ge) ')']);
    xlabel('y(V->D)');ylabel('mean [Dpp]');
end

ScreenSize=get(0,'ScreenSize');
set(gcf,'Position',[ScreenSize(3)*0.02 ScreenSize(4)*0.1 ScreenSize(3)*0.3 ScreenSize(4)*0.4]);
end
function draw_fianl_record(p,re)

subplot(2,2,1);plot(0:p.range_d:p.range_n*p.range_d,re(1:p.range_n+1,:));
title('cytoneme d distribution');
xlabel('cytoneme distance(um)');ylabel('# cytonemes');
for n=1:size(re,2)
    eval(['A3_legend{n}=''' 'g' num2str((n-1)*p.sta_save_period) ''';']);
end
legend(A3_legend);
subplot(2,2,2);plot(re(p.range_n+2,:));
title('cytoneme events during t');
%xlabel('generation');
eval(['xlabel(''' num2str(p.sta_save_period) ' x generation'');']);
ylabel('# of cytoneme events');
subplot(2,2,3);plot(re(p.range_n+3,:));
title('number of cytonemes during t');
eval(['xlabel(''' num2str(p.sta_save_period) ' x generation'');']);
ylabel('# of cytonemes');

subplot(2,2,4);plot(re(p.range_n+4,:));
title('averaged cytoneme distance during t');
eval(['xlabel(''' num2str(p.sta_save_period) ' x generation'');']);
ylabel('averaged cytoneme distance');

ScreenSize=get(0,'ScreenSize');
set(gcf,'Position',[ScreenSize(3)*0.43 ScreenSize(4)*0.48 ScreenSize(3)*0.55 ScreenSize(4)*0.4]);
end
function add_variable(varargin)
if exist('cytoneme.mat')==0
else
load('cytoneme.mat');
end
n=nargin;
for i=1:n
s{i} = inputname(i);
eval([s{i} '=varargin{' num2str(i) '};']);
end
clear s varargin i n
save('cytoneme.mat');

end
function sta=check_cyto(sta,p)

sta.cyto=sortrows(sta.cyto',2)';
sta2.cyto=[];
for j=1:size(sta.cyto,2)
if sta.cyto(2,j)<=p.fig_bound||sta.cyto(6,j)<=p.fig_bound||sta.cyto(2,j)>p.nx-p.fig_bound||sta.cyto(6,j)>p.nx-p.fig_bound
else
    A=sta.cyto(:,j);A(2,1)=A(2,1)-p.fig_bound;A(6,1)=A(6,1)-p.fig_bound;
    sta2.cyto=[sta2.cyto A];
end
end
A=sta2.cyto;sta.cyto=sta2.cyto;
[n,m]=size(A);sta_no=0;
B=zeros(m,p.nx-p.fig_bound*2);
m_inverse=0;
for j=1:m
    if A(6,j)<A(2,j)&&sta_no==0
    m_inverse=j;sta_no=1;
    end
end
m_inverse=0;
if m_inverse==0
    AA=A;
elseif m_inverse==1
    AA=A(:,end:-1:1);
else
    AA=[A(:,1:m_inverse-1) A(:,end:-1:m_inverse)];
end
for j=1:m
    if AA(2,j)<=AA(6,j)
        B(j,AA(2,j))=1; 
        B(j,AA(2,j)+1:AA(6,j))=2; 
    else
        B(j,AA(2,j))=1; 
        B(j,AA(6,j):AA(2,j)-1)=2;     
    end
end
sta.cyto_matrix=B;
end


