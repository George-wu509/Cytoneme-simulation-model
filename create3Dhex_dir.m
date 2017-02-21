function create3Dhex_dir(key,p_pre)
%{ 
 ==== Function description =========================
 Create 3D or 2D hexial mesh data


 George at 2015.10.16 

 ==== vision information ===========================
 revision 2015.09.24 
 > change deg rate from const to const*C
 > p_format=p.gradient_P*p.prob_format*agrad
 > ave_grad=mean(c in angle{}) ---> ave_grad=mean(g in angle{}) in line 638
 > length=min(p.gradient_L*max(agrad,1),max(angle_g{thi}(:,1))); inline 712
 > (0924) add PL=5 mechanism, and fix script2()
 > (1016) add two pathfinding mechanisms
=================================================
%}
if exist('p_pre')==1    
    p=parameter_set(p_pre);
else
    p=parameter_set();  
end
switch key
    case 1
        script1(p) % run cytoneme simulation
    case 2
        script2() % print out figure by sta_save
    case 3
        script3(p) % pitch run cytoneme simulation
    case 4
        script4(p) % run time-lag cytoneme simulation
    case 5
        script5(p) % pitch run time-lag cytoneme simulation
    case 6
        script6(p) % run multi-step cytoneme simulation
    case 7
        script7(p) % run multi-step cytoneme simulation with time-lag
end
end % Create 3D or 2D hexial mesh data

%% 1. Script
function p=parameter_set(p_pre)
% 0. script select
p.name='Model01_1';

p.AN=1;  % pathfinding anglefunction: (1. const),(2. random),(3. max averagd gradient)
p.LE=4; % pathfinding length function: (1. const),(2. random), (3. function of av4raged gradient)(4. highest concentration in thita)
p.PF=4;   % probability of formation function: (1. const),(2. PF const/length), (3. function of averaged gradient),(4. d-effect)
p.MA=2;     % morphogen transport amount function: (1. one cytoneme), (2. N cytoneme)

p.amount=1.5;    % transport morphogen through cytoneme [C/sec] (## no data!!)
p.theta=1; % p.theta=1~12(degree: 0~330) if angle fun=1
p.length=20; % length of const cytoneme value if length fun=1 (## Dpp_avg:20um [242], Hh_avg: 27um in P, 13um in A [26])
p.prob_format=0.2; % probability of cytoneme formation. if prob fun =1 (## estimated!)
p.gradient_L=10; % cytoneme length gradient coefficient. if length fun=3, length=p.gradient_L* averaged gradient(## estimated!)
p.gradient_P2=20; % cytoneme probability gradient coefficient. if probability fun=2, p_format=p.gradient_P*agrad; (## estimated!)
p.gradient_P3=10; % cytoneme probability gradient coefficient. if probability fun=3, p_format=p.gradient_P*agrad; (## estimated!)
p.dfun_option=2; % d-effect options. if probability fun=4 (## estimated)
p.amount_dist=[0,5,10,20;10,5,2,1]; % cutoneme number depend on distance: [4 3 2 1] (## estimated from [26])
p.L_lim=0;

% ----------------------------------------------------------------------

p.max_gener=10;    % max generation
p.morphogen_init=1;     % morphogen: 1=dpp, 2=Hh
p.pitch=1;    % script 3 repeat times

if p.morphogen_init==1 % Dpp cytoneme
    p.prod_init=0;    % max initial morphogen concentration [C]  (## Dpp=0, Hh=0)
    p.prod_time=3.98;    % morphogen production m/t [C/sec]          (## 3.98 mol/um/s)
    p.deg_time=2.52*10^(-4);    % morphogen degradationn m/t [C/sec]         (## Dpp 2.52*10-4 sec-1 [20])
    p.cyto_rate=6; % cytoneme grow rate or transport morphogen rate [L/sec] (## 5-7 um/sec for Dpp [190])
elseif p.morphogen_init==2  %Hh cytoneme
    p.prod_init=0;    % max initial morphogen concentration [C]  (## Dpp=0, Hh=0)
    p.prod_time=0.5;    % morphogen production m/t [C/sec]          (## )
    p.deg_time=3.3*10^(-3);    % morphogen degradationn m/t [C/sec]         (## Hh 3.3*10-3 s-1 [201])
    p.cyto_rate=6; % cytoneme grow rate or transport morphogen rate [L/sec] (## 5-7 um/sec for Dpp [190])
end

% 1. mesh parameters
p.nx=40;   % wing disc nx=120(167um) (## actual width: 223(300um)[242])
p.ny=20;   % wing disc ny=50 (83um)  (## actual width: 111(150um)[242]) 
p.r=0.67;  % 0.67um                  (## actual diameter: 0.5-2.2um[242])
p.w_dppcent=8; % Dpp producing center width=9cells  (## Dpp signaling center wide: 8~10 cells width)
p.height=0;  % if h=-1 means 2D

p.sta_save_period=2;     % figures and record every generations
p.area_rx=5;     % search space x radius (## Dpp:40um =2*Dpp_avg [242],Hh:54um in P,28um in A [26])
p.area_ry=0;     % search space y radius (## Dpp:40um =2*Dpp_avg [242],Hh:54um in P,28um in A [26])
%p.grad_r=40;     % calculate morphogen gradient range radius (## 40um =2*Dpp_avg [242],Hh:54um in P,28um in A [26])

% 3. figure
p.img=1;              % no run any figures sub function
p.save_every_img=0;   % output figures for each display step. 1: output figures, 0: no output figures
p.show_every_img=1;   % show figures for each display step. 1: output figures, 0: no output figures
%p.col_range=[1,max([ceil(p.prod_init),1])];    % jet colormap range
p.col_range=[1,fix(p.prod_time*p.max_gener)];    % jet colormap range
p.alpha1=0.8;    % alpha value for figure1
p.alpha2=0.7;    % alpha value for figure2
p.range_n=5;     % cytoneme distance table value number 
p.range_d=5;     % % cytoneme distance table diff distance

if exist('p_pre')==1    
n=length(p_pre);
for nn=1:n
eval(p_pre{nn});
end
end
end
function script1(p)
tic;
close all;
% set up mesh,init state
disp(['===' p.name '===']);
disp('');
p=create_mesh(p);
sta=initial_state(p);
sta.cyto=zeros(11,1);n=1;re=[];p.i=0;
if p.img==1;
    [sta,p]=show_figure(sta,p,0);
end
sta_save=[];
%sta.cyto_y=zeros(1,p.nx);
% updata generation
for i=1:p.max_gener
    sta=update_apical_dpp(sta,p);p.i=i; % transport update calculation
    sta=avg_y_cyto(sta,p);
    sta_all{i}=sta;
    if rem(i,p.sta_save_period)==0||i==p.max_gener
        disp(['generation  ' num2str(i)]);
        if p.img==1;
        [sta,p]=show_figure(sta,p,i,sta_save);
        end
        sta_save{n}=sta;n=n+1;
        [re,re_txt]=record(sta,p,re);
    end
end

% show Figure 5 during time
if p.img==1;
    F5=figure(5);draw_fianl_record(p,re);
end
save('cytoneme.mat','sta_save','p','re','re_txt','sta_all');
toc;
end % run cytoneme simulation
function script2()
% parameter:
save_gen=0;     % save generation: if value<1 means generation/total, -1=final, 0=all
% ======================
close all
load cytoneme.mat;n=size(sta_save,2);
if save_gen==0
    for i=1:n
        sta=sta_save{i};
        ge=i*p.sta_save_period;
        out_fig(sta,p,re,ge,sta_save);
    end
elseif save_gen==-1
    sta=sta_save{n};
    ge=n*p.sta_save_period;
    out_fig(sta,p,re,ge,sta_save);
else
    sta=sta_save{save_gen};
    ge=save_gen*p.sta_save_period;
    out_fig(sta,p,re,ge,sta_save);
end
T=clock;save([p.name '_result_' num2str(T(2)) '_' num2str(T(3)) '_' num2str(T(4)) '_' num2str(T(5)) '.mat'],'sta_save','p','re','re_txt');
%delete('cytoneme.mat');
end % print out figure by sta_save
function script3(p)

close all;
% set up mesh,init state
disp(['===' p.name '===']);
disp('');
p=create_mesh(p);
sta0=initial_state(p);
%sta0.cyto=zeros(11,1);
%[sta,p]=show_figure(sta,p,0);
% updata generation
for j=1:p.pitch
    tic;
    disp(['[ RUN ' num2str(j) ' ]']);
    sta=sta0;re_pit=[];n=1;
for i=1:p.max_gener
    sta=update_apical_dpp(sta,p); % transport update calculation
    sta=avg_y_cyto(sta,p);
    sta_all{i}=sta;
    if rem(i,p.sta_save_period)==0||i==p.max_gener
        disp(['generation  ' num2str(i)]);
%        [sta,p]=show_figure(sta,p);
        sta_pitch{1,n}.dpp(:,:,j)=sta.dpp;%sta_pitch{1,n}.rgb(:,:,j)=sta.rgb;
        sta_pitch{1,n}.cyto=sta.cyto;
        n=n+1;
        [re_pit,re_txt]=record(sta,p,re_pit);
    end
end
re_pitch(:,:,j)=re_pit;
toc;
end
%sta_save{1,1}=sta0;
for i=1:n-1    
sta_save{1,i}.dpp=mean(sta_pitch{1,i}.dpp,3);
%sta_save{1,i}.rgb=mean(sta_pitch{1,i}.rgb,3);
sta_save{1,i}.cyto=sta.cyto;
end
re=mean(re_pitch,3);
% show Figure 5 during time
save('cytoneme.mat','sta_save','p','re','re_txt','sta_pitch');
if p.img==1;
F5=figure(5);draw_fianl_record(p,re);
end
script2();

end % run cytoneme simulation pitch
function script4(p)
tic;
close all;
% set up mesh,init state
disp(['===' p.name '===']);
disp('');
p=create_mesh(p);
sta=initial_state(p);
sta.cyto=zeros(11,1);n=1;re=[];p.i=0;
if p.img==1;
    [sta,p]=show_figure(sta,p,0);
end
sta_save=[];
%sta.cyto_y=zeros(1,p.nx);
% updata generation
Tlag_list=[];
for i=1:p.max_gener
    [sta,Tlag_list]=update_apical_dpp_time(sta,p,Tlag_list);p.i=i; % transport update calculation
    %sta=update_apical_dpp(sta,p);p.i=i; % transport update calculation
    sta=avg_y_cyto(sta,p);
    sta_all{i}=sta;
    if rem(i,p.sta_save_period)==0||i==p.max_gener
        disp(['generation  ' num2str(i)]);
        if p.img==1;
        [sta,p]=show_figure(sta,p,i,sta_save);
        end
        sta_save{n}=sta;n=n+1;
        [re,re_txt]=record(sta,p,re);
    end
end

% show Figure 5 during time
if p.img==1;
    F5=figure(5);draw_fianl_record(p,re);
end
save('cytoneme.mat','sta_save','p','re','re_txt','sta_all');
toc;
end % run cytoneme simulation TIME LAG
function script5(p)

close all;
% set up mesh,init state
disp(['===' p.name '===']);
disp('');
p=create_mesh(p);
sta0=initial_state(p);
%sta0.cyto=zeros(11,1);
%[sta,p]=show_figure(sta,p,0);
% updata generation
for j=1:p.pitch
    tic;
    disp(['[ RUN ' num2str(j) ' ]']);
    sta=sta0;re_pit=[];n=1;
    Tlag_list=[];
for i=1:p.max_gener
    [sta,Tlag_list]=update_apical_dpp_time(sta,p,Tlag_list);p.i=i; % transport update calculation
    sta=avg_y_cyto(sta,p);
    sta_all{i}=sta;
    if rem(i,p.sta_save_period)==0||i==p.max_gener
        disp(['generation  ' num2str(i)]);
%        [sta,p]=show_figure(sta,p);
        sta_pitch{1,n}.dpp(:,:,j)=sta.dpp;%sta_pitch{1,n}.rgb(:,:,j)=sta.rgb;
        sta_pitch{1,n}.cyto=sta.cyto;
        n=n+1;
        [re_pit,re_txt]=record(sta,p,re_pit);
    end
end
re_pitch(:,:,j)=re_pit;
toc;
end
%sta_save{1,1}=sta0;
for i=1:n-1    
sta_save{1,i}.dpp=mean(sta_pitch{1,i}.dpp,3);
%sta_save{1,i}.rgb=mean(sta_pitch{1,i}.rgb,3);
sta_save{1,i}.cyto=sta.cyto;
end
re=mean(re_pitch,3);
% show Figure 5 during time
save('cytoneme.mat','sta_save','p','re','re_txt','sta_pitch');
if p.img==1;
F5=figure(5);draw_fianl_record(p,re);
end
script2();

end % run cytoneme simulation pitch TIME LAG
function script6(p)
tic;
close all;
% set up mesh,init state
disp(['===' p.name '===']);
disp('');
p=create_mesh(p);
sta=initial_state(p);n=1;re=[];p.i=0;
sta.cyto=[]; %zeros(11,1);
if p.img==1;
    [sta,p]=show_figure(sta,p,0);
end
sta_save=[];Mstep_list=[];Mstep_record=[];
%sta.cyto_y=zeros(1,p.nx);
% updata generation
for i=1:p.max_gener
    [sta,Mstep_list,Mstep_record]=update_apical_dpp_series(sta,p,Mstep_list,Mstep_record);p.i=i; % transport update calculation
    sta=avg_y_cyto(sta,p);
    sta_all{i}=sta;
    if rem(i,p.sta_save_period)==0||i==p.max_gener
        disp(['generation  ' num2str(i)]);
        if p.img==1;
        [sta,p]=show_figure(sta,p,i,sta_save);
        end
        sta_save{n}=sta;n=n+1;
        [re,re_txt]=record(sta,p,re);
    end
end

% show Figure 5 during time
if p.img==1;
    F5=figure(5);draw_fianl_record(p,re);
end
save('cytoneme.mat','sta_save','p','re','re_txt','sta_all');
toc;
end % run cytoneme simulation
function script7(p)
tic;
close all;
% set up mesh,init state
disp(['===' p.name '===']);
disp('');
p=create_mesh(p);
sta=initial_state(p);n=1;re=[];p.i=0;
sta.cyto=[]; %zeros(11,1);
if p.img==1;
    [sta,p]=show_figure(sta,p,0);
end
sta_save=[];Mstep_list=[];Mstep_record=[];
%sta.cyto_y=zeros(1,p.nx);
% updata generation
for i=1:p.max_gener
    [sta,Mstep_list,Mstep_record]=update_apical_dpp_series_time(sta,p,Mstep_list,Mstep_record);p.i=i; % transport update calculation
    sta=avg_y_cyto(sta,p);
    sta_all{i}=sta;
    if rem(i,p.sta_save_period)==0||i==p.max_gener
        disp(['generation  ' num2str(i)]);
        if p.img==1;
        [sta,p]=show_figure(sta,p,i,sta_save);
        end
        sta_save{n}=sta;n=n+1;
        [re,re_txt]=record(sta,p,re);
    end
end

% show Figure 5 during time
if p.img==1;
    F5=figure(5);draw_fianl_record(p,re);
end
save('cytoneme.mat','sta_save','p','re','re_txt','sta_all');
toc;
end % run cytoneme simulation

%% 2. Main program
%2.1 Main process
function p=create_mesh(p)

% create hex mesh parameter
[p.xindex,p.yindex]=meshgrid(1:p.nx,1:p.ny);
p.yc=repmat((1.5*p.r*(p.ny-1):-1.5*p.r:0)',1,p.nx);
p.xc=zeros(p.ny,p.nx);
for wci=1:p.ny
    if rem(wci,2)~=0
        p.xc(p.ny-wci+1,:)=0:2*cosd(30)*p.r:2*cosd(30)*p.r*(p.nx-1);
    else
        p.xc(p.ny-wci+1,:)=-cosd(30)*p.r:2*cosd(30)*p.r:cosd(30)*p.r*(2*p.nx-3);
    end
end

% 2D and 3D hex mesh structure handles
p.h2d_on=0;p.h2dh_on=0;
p.h3d_on=0;p.h3dh_on=0;
p.h2d=zeros(1,p.ny,p.nx);p.h2dh=zeros(1,p.ny,p.nx);
p.h3d=zeros(7,p.ny,p.nx);p.h3dh=zeros(7,p.ny,p.nx);
end % create mesh center location matrix and color map
function sta=initial_state(p)

[ny,nx]=size(p.xindex);
switch p.morphogen_init
    
% (1) dpp concentration on drosophila wing disc
    case 1
sta.dpp=zeros(ny,nx);
sta.dpp(:,fix(nx/2-p.w_dppcent):ceil(nx/2))=p.prod_init;
sta.cyto=[];

% (2) Hh concentration on drosophila wing disc
    case 2
sta.dpp=zeros(ny,nx);
%sta.dpp(:,fix(nx/2-0.5):end)=p.prod_init; 
sta.dpp(:,fix(nx/2+0.5):end)=p.prod_init;   % <-----change - to +
sta.cyto=[];

% (3) output color RGB to new XY_color:EXAMPLE
%{
for co_y=1:ny
    for co_x=1:nx
        XY_color_new.rgb(co_y,co_x)=XY_color_old.xindex(co_y,co_x)+XY_color_old.yindex(co_y,co_x)*(co_y-1);
        if XY_color_new.rgb(co_y,co_x)>size(colormap_now,1);
            XY_color_new.rgb(co_y,co_x)=size(colormap_now,1);
        end
    end
end
%}
end

end % create the initial state
function sta_new=update_apical_dpp(sta_old,p)
sta_new=dpp_grow_reaction(sta_old,p); % morphogen growing and reaction update in every generation
%sta_new.cyto=[];
sta_new.cyto=[];sta_new.cyto_y=[];sta_new.cyto_theta_mean=[];sta_new.cyto_theta_std=[];
sta_new.cyto_length_mean=[];sta_new.cyto_length_std=[];sta_new.dpp_y=[];
% update
    for j=1:p.nx  % j=1:p.nx
    %for j=p.nx:-1:1
        for i=1:p.ny
            if p.xc(i,j)>=0&&p.yc(i,j)>=0
                [i_t,j_t,d,p,p_format]=target_position(sta_old,p,i,j);
                if rand(1,1)<p_format
                    amount=morphogen_amount(d,p);
                    % transport from source to receiving cell
                    if p.morphogen_init==2
                        if j<p.nx/2
                        real_target_dpp=max(sta_old.dpp(i_t,j_t)-amount,0);
                        sta_new.dpp(i_t,j_t)=real_target_dpp;
                        sta_new.dpp(i,j)=sta_old.dpp(i,j)+(sta_old.dpp(i_t,j_t)-real_target_dpp);
                        else
                        real_target_dpp=max(sta_old.dpp(i,j)-amount,0);
                        sta_new.dpp(i,j)=real_target_dpp;
                        sta_new.dpp(i_t,j_t)=sta_old.dpp(i_t,j_t)+(sta_old.dpp(i,j)-real_target_dpp);
                        end
                    else
                        real_target_dpp=max(sta_old.dpp(i_t,j_t)-amount,0);
                        sta_new.dpp(i_t,j_t)=real_target_dpp;
                        sta_new.dpp(i,j)=sta_old.dpp(i,j)+(sta_old.dpp(i_t,j_t)-real_target_dpp);
                    end
                    % update cytoneme events
                    if sta_new.dpp(i_t,j_t)==sta_old.dpp(i_t,j_t)
                    else
                    source_xy=center_coord(p,j,i);target_xy=center_coord(p,j_t,i_t);
                    sta_new.cyto=[sta_new.cyto [i;j;source_xy';i_t;j_t;target_xy';d;amount/p.amount;amount]];
                    end
                    sta_old.dpp=sta_new.dpp;
                end
            end
        end
    end
    for j=1:p.nx
        for i=1:p.ny
            sta_new.dpp(i,j)=sta_new.dpp(i,j)-p.deg_time*sta_new.dpp(i,j);
        end
    end
end 
function [sta_new,Tlag_list]=update_apical_dpp_time(sta_old,p,Tlag_list)
sta_new=dpp_grow_reaction(sta_old,p); % morphogen growing and reaction update in every generation
sta_old=sta_new;
%sta_new.cyto=[];
sta_new.cyto=[];sta_new.cyto_y=[];sta_new.cyto_theta_mean=[];sta_new.cyto_theta_std=[];
sta_new.cyto_length_mean=[];sta_new.cyto_length_std=[];sta_new.dpp_y=[];
% update
    for j=1:p.nx  % j=1:p.nx
    %for j=p.nx:-1:1
        for i=1:p.ny
            if p.xc(i,j)>=0&&p.yc(i,j)>=0
                [i_t,j_t,d,p,p_format]=target_position(sta_old,p,i,j);
                if rand(1,1)<p_format
                    amount=morphogen_amount(d,p);
                    % transport from source to receiving cell
                    if p.morphogen_init==2 %Hh
                        if j<p.nx/2
                        real_target_dpp=max(sta_old.dpp(i_t,j_t)-amount,0);
                        sta_new.dpp(i_t,j_t)=real_target_dpp;
                        sta_new.dpp(i,j)=sta_old.dpp(i,j)+(sta_old.dpp(i_t,j_t)-real_target_dpp);
                        else
                        real_target_dpp=max(sta_old.dpp(i,j)-amount,0);
                        sta_new.dpp(i,j)=real_target_dpp;
                        sta_new.dpp(i_t,j_t)=sta_old.dpp(i_t,j_t)+(sta_old.dpp(i,j)-real_target_dpp);
                        end
                    else %Dpp
                        %real_target_dpp=max(sta_old.dpp(i_t,j_t)-amount,0);
                        %sta_new.dpp(i_t,j_t)=real_target_dpp;
                        %sta_new.dpp(i,j)=sta_old.dpp(i,j)+(sta_old.dpp(i_t,j_t)-real_target_dpp);
                        [Tlag_list]=update_time_list(p,Tlag_list,i,j,i_t,j_t,d,amount);
                    end
                    %update cytoneme events
                    %if sta_new.dpp(i_t,j_t)==sta_old.dpp(i_t,j_t)
                    %else
                    source_xy=center_coord(p,j,i);target_xy=center_coord(p,j_t,i_t);
                    sta_new.cyto=[sta_new.cyto [i;j;source_xy';i_t;j_t;target_xy';d;amount/p.amount;amount]];
                    %end
                    sta_old.dpp=sta_new.dpp;
                end
            end
        end
    end
    [sta_new,Tlag_list]=update_time_transport(p,Tlag_list,sta_new);
    for j=1:p.nx
        for i=1:p.ny
            sta_new.dpp(i,j)=max(sta_new.dpp(i,j)-p.deg_time,0);
        end
    end
end % transport update calculation
function [sta_new,Mstep_list,Mstep_record]=update_apical_dpp_series(sta_old,p,Mstep_list,Mstep_record)
sta_new=dpp_grow_reaction(sta_old,p); % morphogen growing and reaction update in every generation
%sta_new.cyto=[];
sta_new.cyto=[];sta_new.cyto_y=[];sta_new.cyto_theta_mean=[];sta_new.cyto_theta_std=[];
sta_new.cyto_length_mean=[];sta_new.cyto_length_std=[];sta_new.dpp_y=[];c_i=1;
% update
    for j=1:p.nx  % j=1:p.nx
    %for j=p.nx:-1:1
        for i=1:p.ny
            if p.xc(i,j)>=0&&p.yc(i,j)>=0
                [i_t,j_t,d,p,p_format]=target_position(sta_old,p,i,j);
                if rand(1,1)<p_format
                    [Mstep_list,Mstep_record,c_i]=multi_transport_record(Mstep_list,Mstep_record,i,j,i_t,j_t,c_i,d);
                else
                    [Mstep_list,Mstep_record,sta_new]=multi_transport_update(Mstep_list,Mstep_record,sta_new,i,j,p);
                end
            end
        end
    end
    for j=1:p.nx
        for i=1:p.ny
            sta_new.dpp(i,j)=max(sta_new.dpp(i,j)-p.deg_time,0);
        end
    end
end 
function [sta_new,Mstep_list,Mstep_record]=update_apical_dpp_series_time(sta_old,p,Mstep_list,Mstep_record)
sta_new=dpp_grow_reaction(sta_old,p); % morphogen growing and reaction update in every generation
%sta_new.cyto=[];
sta_new.cyto=[];sta_new.cyto_y=[];sta_new.cyto_theta_mean=[];sta_new.cyto_theta_std=[];
sta_new.cyto_length_mean=[];sta_new.cyto_length_std=[];sta_new.dpp_y=[];c_i=1;Tlag_list=[];
% update
    for j=1:p.nx  % j=1:p.nx
    %for j=p.nx:-1:1
        for i=1:p.ny
            if p.xc(i,j)>=0&&p.yc(i,j)>=0
                [i_t,j_t,d,p,p_format]=target_position(sta_old,p,i,j);
                if rand(1,1)<p_format
                    [Mstep_list,Mstep_record,c_i]=multi_transport_record(Mstep_list,Mstep_record,i,j,i_t,j_t,c_i,d);
                else
                    [Mstep_list,Mstep_record,sta_new,amount]=multi_transport_update_time(Mstep_list,Mstep_record,sta_new,i,j,p);
                    if amount~=0
                        [Tlag_list]=update_time_list_series(p,Tlag_list,i,j,i_t,j_t,d,amount);
                    end
                end
            end
        end
    end
    for j=1:p.nx
        for i=1:p.ny
            sta_new.dpp(i,j)=max(sta_new.dpp(i,j)-p.deg_time,0);
        end
    end
    [sta_new,Tlag_list]=update_time_transport(p,Tlag_list,sta_new);
end 
function sta_new=dpp_grow_reaction(sta_old,p)
[ny,nx]=size(sta_old.dpp);
sta_new=sta_old;
% For growing part
switch p.morphogen_init   
% (1) dpp concentration on drosophila wing disc
case 1
sta_grow=zeros(ny,nx);
sta_grow(:,fix(nx/2-p.w_dppcent):ceil(nx/2))=p.prod_time;
sta_new.dpp=sta_new.dpp+sta_grow;
% (2) Hh concentration on drosophila wing disc
case 2
sta_grow=zeros(ny,nx);
sta_grow(:,fix(nx/2-0.5):end)=p.prod_time;
sta_new.dpp=sta_new.dpp+sta_grow;
end

% For reaction part
sta_reaction=zeros(ny,nx);
sta_new.dpp=sta_new.dpp+sta_reaction;
end % morphogen growing and reaction update

%% 3. Sub function
%3.1 cytoneme transport function
function [i_t,j_t,d,p,p_format]=target_position(sta_old,p,i,j)
% [i_t,j_t,d,p,dmr,dms]=target_position(sta_old,p,i,j)
% 20140513: dont consider the dist effect on TP(), so dont output deff and meff
    % define search sapce for receiving cell
    i_o=max(i-p.area_ry,1);j_o=max(j-p.area_rx,1);
    se_regi=sta_old.dpp(i_o:min(i+p.area_ry,p.ny),j_o:min(j+p.area_rx,p.nx));
    [angle_g,ave_grad]=direction_gradient(se_regi,p,i,j,i_o,j_o); % 20140817 add direction_gradient() here
    [i_t,j_t,d,p_format]=target_find(angle_g,ave_grad,sta_old,p,i,j); % 20140817 add direction_gradient() here
end % TP() - target morphogen source cell()
function [angle_g,ave_grad]=direction_gradient(se_regi,p,i,j,i_o,j_o)
[a,b]=size(se_regi);
angle_table=zeros(a*b,7); % angle_table=[thi d i j x y c]
x_orig=p.xc(i,j);y_orig=p.yc(i,j);
for ii=1:a  % only ii,jj are se_regi coordination system
    for jj=1:b
        x_target=p.xc(ii+i_o-1,jj+j_o-1);y_target=p.yc(ii+i_o-1,jj+j_o-1);
        if x_target-x_orig<0
            angle_table(jj+(ii-1)*b,1)=atand((y_target-y_orig)/(x_target-x_orig))+180;
        elseif y_target-y_orig<0
            angle_table(jj+(ii-1)*b,1)=atand((y_target-y_orig)/(x_target-x_orig))+360;
        else
            if abs(y_target-y_orig)<0.001&&abs(x_target-x_orig)<0.001
                angle_table(jj+(ii-1)*b,1)=0;
            else          
                angle_table(jj+(ii-1)*b,1)=atand((y_target-y_orig)/(x_target-x_orig));
            end
        end
        angle_table(jj+(ii-1)*b,2)= distance(x_orig,y_orig,x_target,y_target); %dist between center to target
        angle_table(jj+(ii-1)*b,3)=ii+i_o-1; % ii(search space) --> i(target in disc plane)
        angle_table(jj+(ii-1)*b,4)=jj+j_o-1; % jj(search space) --> j(target in disc plane)
        angle_table(jj+(ii-1)*b,5)=x_target;angle_table(jj+(ii-1)*b,6)=y_target; % real coordinates for target
        angle_table(jj+(ii-1)*b,7)=se_regi(ii,jj); % morphogen concentration
    end
end
angle_table=sortrows(angle_table); % all elemetns in search space sorted by angle
ave_grad=zeros(1,12);
q=1; % thita= 0
% angle_t: temp [thi d i j x y c]
%angle_t0=[0 0 i j x_orig y_orig se_regi(i-i_o+1,j-j_o+1)];
c0=se_regi(i-i_o+1,j-j_o+1);

angle_t=angle_table(angle_table(:,1)>=345|angle_table(:,1)<15,:); % angle_t: part of 'angle_table' for range of angle
angle_t=angle_t(angle_t(:,2)~=0,:);
angle_t=sortrows(angle_t,2);
angle_g{q}=position_to_gradient(angle_t,q,c0); % angle_g: re-sorted and ranked from 'angle_t' using distance
ave_grad(1,1)=mean(angle_g{1}(:,3)); % <------old ave gradient(07/08)
%ave_grad(1,1)=mean(angle_g{1}(:,2)); % <------new ave gradient(07/08)
for q=2:12 % thita=30,60...330
    angle_t=angle_table(angle_table(:,1)>=(q-1)*30-15&angle_table(:,1)<(q-1)*30+15,:);
    if isempty(angle_t)==1
        angle_g{q}=[];
        ave_grad(1,q)=0;
    else
        angle_t=sortrows(angle_t,2);
        angle_g{q}=position_to_gradient(angle_t,q,c0);
        ave_grad(1,q)=mean(angle_g{q}(:,3)); % <------old ave gradient(07/08)
        %ave_grad(1,q)=mean(angle_g{q}(:,2)); % <------new ave gradient(07/08)
    end
end
end
function [i_t,j_t,d,p_format]=target_find(angle_g,ave_grad,sta_old,p,i,j)
stop=0;
% pathfinding angle function - thi:1~12, ag= averaged gradient of thi
if p.AN==1 % pathfinding angle (1. const)
    if j<(p.nx-p.w_dppcent)/2
        thi=p.theta;
    else
        thi=p.theta+6;
    end
    ag=ave_grad(thi);
elseif p.AN==2 % pathfinding angle (2. random)
    thi=min(max(ceil(12*rand(1)),0),size(ave_grad,2));
    ag=ave_grad(thi);
elseif p.AN==3 % pathfinding angle (3. max averaged gradient)
    if p.morphogen_init==1
        [ag,thi]=max(ave_grad);
    elseif p.morphogen_init==2
        if j<p.nx/2
          [ag,thi]=max(ave_grad);
        else
          [ag,thi]=min(ave_grad);  
        end
    end
end
if ag==0
    stop=1;
    agrad=0;
    length_maxc=0;
else
    agrad=ave_grad(1,thi); % agrad: averaged gradient
    if p.morphogen_init==1
        [~,b]=max(angle_g{thi}(:,2));
        length_maxc=angle_g{thi}(b,1);
    elseif p.morphogen_init==2
        if j<p.nx/2
            [~,b]=max(angle_g{thi}(:,2));
            length_maxc=angle_g{thi}(b,1);
        else
            [~,b]=min(angle_g{thi}(:,2));
            length_maxc=angle_g{thi}(b,1);
        end
    end
    if agrad<0
        agrad=-agrad;
    end
end

if isempty(angle_g{thi})==1
    d=0;stop=1;
    i_t=i;j_t=j;
else
    % pathfinding length function
    if p.LE==1&&p.L_lim==1          % pathfinding length (1, const)
        length=min(max(ceil(p.length),1),max(angle_g{thi}(:,1)));
    elseif p.LE==2&&p.L_lim==1        % pathfinding length (2, rand)
        length=min(max(ceil(p.length*rand(1)),1),max(angle_g{thi}(:,1)));
    elseif p.LE==3&&p.L_lim==1        % pathfinding length (3, fun of aved gradient)
        length=min(p.gradient_L*max(agrad,1),max(angle_g{thi}(:,1)));
    elseif p.LE==4&&p.L_lim==1        % pathfinding length (4, highest concentration in thita)
        length=length_maxc;
    elseif p.LE==5&&p.L_lim==1        % pathfinding length (5, fun of inverse gradient)
        length=min(p.gradient_L/max(agrad,1),max(angle_g{thi}(:,1)));
    elseif p.LE==1&&p.L_lim==0        % pathfinding length (1, const)no limit
        length=max(ceil(p.length),1);
    elseif p.LE==2&&p.L_lim==0        % pathfinding length (2, rand)no limit
        length=max(ceil(p.length*rand(1)),1);
    elseif p.LE==3&&p.L_lim==0        % pathfinding length (3, fun of aved gradient)no limit
        length=p.gradient_L*max(agrad,1);
    elseif p.LE==4&&p.L_lim==0        % pathfinding length (4, highest concentration in thita)no limit
        length=length_maxc;
    elseif p.LE==5&&p.L_lim==0       % pathfinding length (5, fun of inverse gradient)no limit
        length=p.gradient_L/max(agrad,1);
    end
    angle_gg=angle_g{1,thi}(angle_g{1,thi}(:,1)>0,:);
    [~,n]=min(abs(angle_gg(:,1)-length));
    i_t=angle_gg(n,4);j_t=angle_gg(n,5);
    d=angle_gg(n,1);
end
if p.morphogen_init==1
    if sta_old.dpp(i_t,j_t)==0
        stop=1;
    else
        stop=0;
    end
end
% probability of formation function 
if p.PF==1&&stop~=1 % probability of formation (1. const)
    p_format=p.prob_format;
elseif p.PF==2&&stop~=1 % probability of formation (2. PF const/length)
    p_format=p.gradient_P2*probability_formation(thi,d,p);
elseif p.PF==3&&stop~=1 % probability of formation (3. fun of aved gradient)
    p_format=p.gradient_P3*p.prob_format*agrad;
elseif p.PF==4&&stop~=1 % probability of formation (3. p=d_effect)
    p_format=dist_effect(d,p);
elseif p.PF==5&&stop~=1 % probability of formation (3. p=d_effect)
    p_format=p.gradient_P3*p.prob_format/agrad;
else
    p_format=0;
end
end
function amount=morphogen_amount(d,p)
switch p.MA
    case 1 % one event one cytoneme, amount=M
        amount=p.amount;
    case 2 % amount=NCxM
        nn=max(find(d>=p.amount_dist(1,:)));
        amount=p.amount_dist(2,nn)*p.amount; % decide the cytoneme transport amount
    case 3 % time delay(cytoneme grow, transport)
end
end % MA() - morphogen amount for cytoneme transport()

function angle_g=position_to_gradient(angle_t,q,c0)
nn=ceil(max(angle_t(:,2)));
angle_g=zeros(nn,5);
if q==1  % thita=0 tranfrom 345~360degree to 0~15
    q0=size(angle_t,1);
    for qq0=1:q0
       if angle_t(qq0,1)>345
           angle_t(qq0,1)=360-angle_t(qq0,1);
       end
    end
end
for n=1:nn % in specific thi, from distance=1 to nn
    temp=angle_t(angle_t(:,2)>=n-1&angle_t(:,2)<n,:); %distance from n-1 to n
    [~,b]=min(abs(temp(:,1)-(q-1)*30));c=mean(temp(:,7));
    if n==1
        g=c-c0; % g= c(n)-c(n-1), and distance(n,n-1)=1
    else
        g=c-angle_g(n-1,2);
    end
    % g=(angle_t(n,7)-angle_t(n-1,7))/(distance(angle_t(n-1,5),angle_t(n-1,6),angle_t(n,5),angle_t(n,6)));
    if isempty(temp)
        angle_g(n,:)=zeros(1,5);
    else    
        angle_g(n,:)=[n c g temp(b,3) temp(b,4)];
    end
end
end
function d=grad_to_length(a)
d=5*a;
end % averaged gradient-cytoneme length function
function p_format=probability_formation(~,length,p)
p_format=p.prob_format*(1/length);
end % formation probability function
function d_effect=dist_effect(d,p)
% fist effect function parameters
p.center=20;     % p_thread normpdf(d,p.center,p.the)*p.thread_multi for case2
p.the=10;     % p_thread normpdf(d,p.center,p.the)*p.thread_multi
p.d_coef=0.5;
% -------------------------------------------------
switch p.dfun_option
    case 1
        d_effect = p.d_coef;
    case 2
        xmax = normpdf(p.center,p.center,p.the);
        d_effect = p.d_coef*normpdf(d,p.center,p.the)/xmax*2;
    case 3
        d_effect = p.d_coef*lognpdf(d,log(25),0.5)/0.0362*2;
end
end
function m_effect=mor_effect(mr,ms,dmr,dms,p,m_dist)
% m_effect= [morphogen] effect on cytoneme formation from 0(lowest)~2(highest)
switch m_dist
    case 1
        m_effect = 1;
    case 2  % delta(MS-MR)
        m_effect = max(2*(ms-mr)/p.prod_init,0);
    case 3  % morphogen gradient effect on morphogen receptor
        m_effect = max(abs(dmr/p.amount),0);
    case 4 % morphogen gradient effect on morphogen source
        m_effect = max(dms/p.amount,0);
end
end
function sta=avg_y_cyto(sta,p)
sta.cyto_y=zeros(1,p.nx);
all_ttt{p.nx}=[];all_ddd{p.nx}=[];

[~,cn]=size(sta.cyto);
for i=1:cn
    c=sta.cyto(1,i);d=sta.cyto(5,i); %c(x_o), d(x_t)
    a=sta.cyto(2,i);b=sta.cyto(6,i); %a(y_o), b(y_t)
    
    if c==d
        if b-a>=0
            ttt=0;
        else
            ttt=180;
        end
    else
        if b-a>=0
            if d-c>=0
                ttt=90-atand((b-a)/(d-c));
            else
                ttt=-(90+atand((b-a)/(d-c)));
            end
        else
            ttt=180+atand((b-a)/(d-c));
        end
    end    
    ddd=sta.cyto(9,i);
      
    if a<=b
        A=a:b;
    else
        A=b:a;
    end
    for j=A
        sta.cyto_y(1,j)=sta.cyto_y(1,j)+1;
        all_ttt{j}=[all_ttt{j} ttt];
        all_ddd{j}=[all_ddd{j} ddd];
    end
end
for k=1:p.nx
    sta.cyto_theta_mean(k)=mean(all_ttt{k});sta.cyto_theta_std(k)=std(all_ttt{k});
    sta.cyto_length_mean(k)=mean(all_ddd{k});sta.cyto_length_std(k)=std(all_ddd{k});
end
sta.dpp_y=mean(sta.dpp);

sta=check_cyto(sta,p);
end
function [Tlag_list]=update_time_list(p,Tlag_list,i,j,i_t,j_t,d,amount)
[t1,t2]=distance_time_fun(d,p);
if p.morphogen_init==1
    Tlag_one=[t1,i,j,t2,amount,i_t,j_t,t1,-amount];
    Tlag_list=[Tlag_list;Tlag_one];
elseif p.morphogen_init==2
    Tlag_one=[t2,i,j,t2,-amount,i_t,j_t,t1,amount];
    Tlag_list=[Tlag_list;Tlag_one];
end

end
function [sta_new,Tlag_list]=update_time_transport(p,Tlag_list,sta_old)
Tlag_list=sortrows(Tlag_list);[a,~]=size(Tlag_list);sta_new=sta_old;
if p.morphogen_init==1 %dpp
    for i=1:a
        if Tlag_list(i,1)==0
            amou=sta_old.dpp(Tlag_list(i,6),Tlag_list(i,7))-max(sta_old.dpp(Tlag_list(i,6),Tlag_list(i,7))+Tlag_list(i,9),0); % amou: real amount of morphogen transport
            Tlag_list(i,5)=amou;Tlag_list(i,9)=-amou;            
        end
        if Tlag_list(i,8)==0
            sta_new.dpp(Tlag_list(i,6),Tlag_list(i,7))=sta_old.dpp(Tlag_list(i,6),Tlag_list(i,7))+Tlag_list(i,9);
        end
        if Tlag_list(i,4)==0
            sta_new.dpp(Tlag_list(i,2),Tlag_list(i,3))=sta_old.dpp(Tlag_list(i,2),Tlag_list(i,3))+Tlag_list(i,5);
        end
        sta_old.dpp=sta_new.dpp;
        if Tlag_list(i,4)<=0&&Tlag_list(i,8)<=0
            Tlag_list(i,1)=-100;           
        end       
        Tlag_list(i,1)=Tlag_list(i,1)-1;Tlag_list(i,4)=Tlag_list(i,4)-1;Tlag_list(i,8)=Tlag_list(i,8)-1;
    end
    if isempty(Tlag_list)==1
    else
    Tlag_list=Tlag_list(Tlag_list(:,1)>-80,:);
    end
else %hh
    for i=1:a
        if Tlag_list(i,1)==0
            amou=sta_old.dpp(Tlag_list(i,2),Tlag_list(i,3))-max(sta_old.dpp(Tlag_list(i,2),Tlag_list(i,3))+Tlag_list(i,5),0); % amou: real amount of morphogen transport
            Tlag_list(i,9)=amou;Tlag_list(i,5)=-amou;            
        end
        if Tlag_list(i,8)==0
            sta_new.dpp(Tlag_list(i,6),Tlag_list(i,7))=sta_old.dpp(Tlag_list(i,6),Tlag_list(i,7))+Tlag_list(i,9);
        end
        if Tlag_list(i,4)==0
            sta_new.dpp(Tlag_list(i,2),Tlag_list(i,3))=sta_old.dpp(Tlag_list(i,2),Tlag_list(i,3))+Tlag_list(i,5);
        end
        sta_old.dpp=sta_new.dpp;
        if Tlag_list(i,4)<=0&&Tlag_list(i,8)<=0
            Tlag_list(i,1)=-100;           
        end       
        Tlag_list(i,1)=Tlag_list(i,1)-1;Tlag_list(i,4)=Tlag_list(i,4)-1;Tlag_list(i,8)=Tlag_list(i,8)-1;
    end
    if isempty(Tlag_list)==1
    else
    Tlag_list=Tlag_list(Tlag_list(:,1)>-80,:);
    end
end

end
function [Tlag_list]=update_time_list_series(p,Tlag_list,i,j,i_t,j_t,d,amount)
[t1,t2]=distance_time_fun(d,p);t1=0;
if p.morphogen_init==1
    Tlag_one=[t1,i,j,t2,amount,i_t,j_t,t1,-amount];
    Tlag_list=[Tlag_list;Tlag_one];
elseif p.morphogen_init==2
    Tlag_one=[t2,i,j,t2,-amount,i_t,j_t,t1,amount];
    Tlag_list=[Tlag_list;Tlag_one];
end

end
function [t1,t2]=distance_time_fun(d,p)
if p.morphogen_init==1
    t1=max(round(d/p.cyto_rate),0);
    t2=max(round(d/p.cyto_rate),0)+t1;
elseif p.morphogen_init==2
    t1=max(round(d/p.cyto_rate),0);
    t2=0;
end
end
function sta=check_cyto(sta,p)

A=sta.cyto;[n,m]=size(A);
B=zeros(m,p.nx);    
for j=1:m
    if A(2,j)<=A(6,j)
        B(j,A(2,j))=1; 
        B(j,A(2,j)+1:A(6,j))=2; 
    else
        B(j,A(2,j))=1; 
        B(j,A(6,j):A(2,j)-1)=2;     
    end
end
sta.cyto_matrix=B;
end
function [Mstep_list,Mstep_record,c_i]=multi_transport_record(Mstep_list,Mstep_record,i,j,i_t,j_t,c_i,d)
if isempty(Mstep_list)==1
   Mstep_list=[c_i;i;j;i_t;j_t;1;d];
   Mstep_record{c_i}=[i j;i_t j_t];
   c_i=c_i+1;
else
m=find(Mstep_list(4,:)==i&Mstep_list(5,:)==j, 1);
if isempty(m)==1
    Mstep_list_temp=[c_i;i;j;i_t;j_t;1;d];
    Mstep_list=[Mstep_list Mstep_list_temp];
    Mstep_record{c_i}=[i j;i_t j_t];
    c_i=c_i+1;
else
    Mstep_list(4,m)=i_t;Mstep_list(5,m)=j_t;
    Mstep_list(6,m)=Mstep_list(6,m)+1;Mstep_list(7,m)=Mstep_list(7,m)+d;
    Mstep_record{Mstep_list(1,m)}=[Mstep_record{Mstep_list(1,m)};[i_t j_t]];
end
end
end
function [Mstep_list,Mstep_record,sta_new]=multi_transport_update(Mstep_list,Mstep_record,sta_old,i,j,p)
sta_new=sta_old;
if isempty(Mstep_list)==1
else
    m=find(Mstep_list(4,:)==i&Mstep_list(5,:)==j, 1);
    if isempty(m)==1
    else
        amount=morphogen_amount(Mstep_list(7,m),p);
        i_o=Mstep_list(2,m);j_o=Mstep_list(3,m);
        % transport from source to receiving cell
        if p.morphogen_init==2
            if j<p.nx/2
            real_target_dpp=max(sta_old.dpp(i,j)-amount,0);
            sta_new.dpp(i,j)=real_target_dpp;
            sta_new.dpp(i_o,j_o)=sta_old.dpp(i_o,j_o)+(sta_old.dpp(i,j)-real_target_dpp);
            else
            real_target_dpp=max(sta_old.dpp(i_o,j_o)-amount,0);
            sta_new.dpp(i_o,j_o)=real_target_dpp;
            sta_new.dpp(i,j)=sta_old.dpp(i,j)+(sta_old.dpp(i_o,j_o)-real_target_dpp);
            end
        else
            real_target_dpp=max(sta_old.dpp(i,j)-amount,0);
            sta_new.dpp(i,j)=real_target_dpp;
            sta_new.dpp(i_o,j_o)=sta_old.dpp(i_o,j_o)+(sta_old.dpp(i,j)-real_target_dpp);
        end                  
        % update cytoneme events
        if sta_new.dpp(i,j)==sta_old.dpp(i,j)
        else
        source_xy=center_coord(p,j_o,i_o);target_xy=center_coord(p,j,i);
        sta_new.cyto=[sta_new.cyto [i_o;j_o;source_xy';i;j;target_xy';Mstep_list(7,m);amount/p.amount;amount]];
        end
        mmm=size(Mstep_list,2);
        if m==1
            Mstep_list=Mstep_list(:,2:end);
        elseif m==mmm
            Mstep_list=Mstep_list(:,1:m-1);
        else
            Mstep_list=[Mstep_list(:,1:m-1) Mstep_list(:,m+1:end)];
        end
        step_record{m}=[];
    end
end
end
function [Mstep_list,Mstep_record,sta_new,amount]=multi_transport_update_time(Mstep_list,Mstep_record,sta_old,i,j,p)
sta_new=sta_old;
if isempty(Mstep_list)==1
    amount=0;
else
    m=find(Mstep_list(4,:)==i&Mstep_list(5,:)==j, 1);
    if isempty(m)==1
        amount=0;
    else
        amount=morphogen_amount(Mstep_list(7,m),p);
        i_o=Mstep_list(2,m);j_o=Mstep_list(3,m);
        % transport from source to receiving cell
        if p.morphogen_init==2
            if j<p.nx/2
            real_target_dpp=max(sta_old.dpp(i,j)-amount,0);
            sta_new.dpp(i,j)=real_target_dpp;
            %sta_new.dpp(i_o,j_o)=sta_old.dpp(i_o,j_o)+(sta_old.dpp(i,j)-real_target_dpp);
            else
            real_target_dpp=max(sta_old.dpp(i_o,j_o)-amount,0);
            %sta_new.dpp(i_o,j_o)=real_target_dpp;
            sta_new.dpp(i,j)=sta_old.dpp(i,j)+(sta_old.dpp(i_o,j_o)-real_target_dpp);
            end
        else
            real_target_dpp=max(sta_old.dpp(i,j)-amount,0);
            sta_new.dpp(i,j)=real_target_dpp;
            real_transport=(sta_old.dpp(i,j)-real_target_dpp);
            %sta_new.dpp(i_o,j_o)=sta_old.dpp(i_o,j_o)+(sta_old.dpp(i,j)-real_target_dpp);
        end
        amount=real_target_dpp;
        % update cytoneme events
        if sta_new.dpp(i,j)==sta_old.dpp(i,j)
        else
        source_xy=center_coord(p,j_o,i_o);target_xy=center_coord(p,j,i);
        sta_new.cyto=[sta_new.cyto [i_o;j_o;source_xy';i;j;target_xy';Mstep_list(7,m);amount/p.amount;amount]];
        end
        mmm=size(Mstep_list,2);
        if m==1
            Mstep_list=Mstep_list(:,2:end);
        elseif m==mmm
            Mstep_list=Mstep_list(:,1:m-1);
        else
            Mstep_list=[Mstep_list(:,1:m-1) Mstep_list(:,m+1:end)];
        end
        step_record{m}=[];
    end
end
end

%3.2 figure formation
function [sta,p]=show_figure(sta,p,ge,sta_save)
    if p.show_every_img==1;
        % draw mesh color distribution
        F1=figure('visible','on');
        [sta,p]=draw_mesh(sta,p,ge);
        % draw 3D mesh color distribution
        F2=figure('visible','on');drawsurf(sta,p,ge);
    else
        % draw mesh color distribution
        F1=figure('visible','off');[sta,p]=draw_mesh(sta,p,ge);
        % draw 3D mesh color distribution
        F2=figure('visible','off');drawsurf(sta,p,ge);
    end
    if p.save_every_img==1;
        eval(['saveas(F1,''f1_' num2str(p.i) ''',''jpg'');']);
        eval(['saveas(F2,''f2_' num2str(p.i) ''',''jpg'');']);
    end
    p=draw_cytoneme(sta,p,ge);
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
end % draw hex mesh figure
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
function out_fig(sta,p,re,ge,sta_save)
    % draw mesh color distribution
    F1=figure('visible','off');p.height=0;sta=draw_mesh(sta,p,ge);
    % draw 3D mesh color distribution
    F2=figure('visible','off');drawsurf(sta,p,ge);
    F3=figure('visible','off');p=draw_cytoneme(sta,p,ge,F3);
    %{
    if exist('sta_save','var')~=1
        F4=figure(4);draw_section(sta,ge);
    else
        F4=figure(4);draw_section(sta,ge,sta_save);
    end
    F5=figure(5);draw_fianl_record(p,re);
    %}
    eval(['saveas(F1,''F1_' num2str(ge) ''', ''fig'');']);
    eval(['saveas(F2,''F2_' num2str(ge) ''', ''fig'');']);
    eval(['saveas(F3,''F3_' num2str(ge) ''', ''fig'');']);
    eval(['print(F1,''-r' num2str(300) ''',''' '-dtiff' ''',''' 'F1_' num2str(ge) '.' 'tif' ''');']);
    eval(['print(F2,''-r' num2str(300) ''',''' '-dtiff' ''',''' 'F2_' num2str(ge) '.' 'tif' ''');']);
    eval(['print(F3,''-r' num2str(300) ''',''' '-dtiff' ''',''' 'F3_' num2str(ge) '.' 'tif' ''');']);
end
function p=draw_cytoneme(sta,p,ge,F3)
%
p.show_every_img=1;
if p.show_every_img==1;
    F3=figure('visible','on');
else
    F3=figure('visible','off');
end
%
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
        eval(['saveas(F3,''f3_' num2str(p.i) ''',''jpg'');']);
    end
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
%3.3 oter function
function xy=center_coord(p,jj,ii)
xy=[p.xc(ii,jj) p.yc(ii,jj)];
end
function [re,re_txt]=record(sta,p,re)
re_txt={'number of cytonemes in distance range groups';'#cytoneme events';'# cytonemes';'mean cytoneme distances'};    
A=zeros(p.range_n+4,1);
if isempty(sta.cyto)==1
else
for j=1:p.range_n
    A(j,1)=size(find(sta.cyto(9,:)>=p.range_d*(j-1)&sta.cyto(9,:)<p.range_d*j),2);
end
A(p.range_n+1,1)=size(find(sta.cyto(9,:)>=p.range_d*j),2);
A(p.range_n+2,1)=size(find(sta.cyto(10,:)),2);A(p.range_n+3,1)=sum(sta.cyto(10,:));
A(p.range_n+4,1)=mean(sta.cyto(9,:),2);
end
if exist('re','var')==0
    re=A;
else
    re=[re A];
end
end
%3.4 backup code
function p_thread=cytoneme_formation(d,p,mr,ms,dmr,dms)
switch p.CF
    case 1 % single constant threshold
        p_thread=p.CF_base;
    case 2 % cytoneme formation dist distribution effect
        p_thread=p.CF_base*dist_effect(d,p,p.d_dist);
    case 3 % delta(MS-MR)
        p_thread=p.CF_base*mor_effect(mr,ms,dmr,dms,p,p.m_dist);
    case 4 % morphogen spatial distribution effect
        p_thread=p.CF_base*dist_effect(d,p,p.d_dist)*mor_effect(mr,ms,dmr,dms,p,p.m_dist);
end
end % CF() - cytoneme formation event()