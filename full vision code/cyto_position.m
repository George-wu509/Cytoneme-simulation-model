function cyto_position()

%----------------------------
generation=1;
T_name=[52:60];
O_name=[52:60];
% ---------------------------
load cytoneme.mat;
Data=sta_all{1,generation}.cyto;

if isempty(T_name)~=1
for i=1:size(T_name,2)    
    T_list_temp=[];
    for k=1:size(Data,2)
        if Data(6,k)==T_name(i)
            T_list_temp=[T_list_temp;[Data(2,k) T_name(i)]];
        end      
    end
    T_list{i}=T_list_temp;
end
end

if isempty(O_name)~=1
for i=1:size(O_name,2)  
    O_list_temp=[];
    for k=1:size(Data,2)
        if Data(2,k)==O_name(i)
            O_list_temp=[O_list_temp;[O_name(i) Data(6,k)]];
        end      
    end
    O_list{i}=O_list_temp;
end
end
if isempty(T_list)~=1
    creatF(T_list,2);
end
if isempty(O_list)~=1
    creatF(O_list,1)
end
save('cyto_loca.mat','T_list','O_list')
end

function creatF(Data,aa)
% Create figure
figure1=figure;

nn=size(Data,2);
Datamax=0;Datamin=200;

for ii=1:nn
    Datamax=max(max(Data{1,ii}(:,3-aa)),Datamax);
    Datamin=min(min(Data{1,ii}(:,3-aa)),Datamin);
end

for ii=1:nn
% Create subplot
eval(['subplot' num2str(ii) ' = subplot(' num2str(nn) ',1,' num2str(ii) ',''Parent'',figure' num2str(1) ',''CLim'',[1 2]);']);
eval(['box(subplot' num2str(ii) ',''on'');']);


% Create patch
%patch('Parent',subplot1,'VertexNormals',VertexNormals1,'YData',YData1,'XData',XData1,'Vertices',Vertices1,'Faces',Faces1,'FaceColor','flat','FaceVertexCData',FaceVertexCData1,'CData',CData1);
%eval(['patch(''Parent'',subplot' num2str(ii) ',''YData'',Data{1,' num2str(ii) '}(:,3-aa),''XData'',Datamin:1:Datamax);']);
eval(['hist(Data{1,' num2str(ii) '}(:,' num2str(3-aa) '),' num2str(Datamin) ':' num2str(Datamax) ')']);
if aa==2
    eval(['ylabel(''T= ' num2str(Data{1,ii}(1,aa)) ''');']);
elseif aa==1
    eval(['ylabel(''O= ' num2str(Data{1,ii}(1,aa)) ''');']);
end



end
end