function fMenuOffsetMap(func,varargin)
switch func
    case 'AddTo'
        AddTo(varargin{1});
    case 'Clear'
        Clear(varargin{1});
    case 'Match'
        Match(varargin{1});
    case 'MatchAutomatic'
        MatchAutomatic(varargin{1});
    case 'Load'
        Load(varargin{1});
    case 'Save'
        Save(varargin{1});    
    case 'Correct'
        Correct(varargin{1});          
    case 'Show'
        Show;          
    case 'AlignCheck'
        AlignCheck;
    case 'CreateOffsetMap'
        CreateOffsetMap;
    case 'Apply'
        ApplyOffset;
end

function AlignCheck
global Molecule;
global Filament;
if strcmp(get(gcbo,'Checked'),'on')==1
    set(gcbo,'Checked','off');
    mode = 1;
else
    set(gcbo,'Checked','on');
    mode = 0;
end
Molecule = fTransformCoord(Molecule,mode,0);
Filament = fTransformCoord(Filament,mode,1);
fShow('Image');
fShow('Tracks');

function CreateOffsetMap
global Molecule;
global Filament;
hMainGui = getappdata(0,'hMainGui');
set(hMainGui.Menu.mAlignChannels,'Checked','off');
Molecule = fTransformCoord(Molecule,1,0);
Filament = fTransformCoord(Filament,1,1);
fShared('UpdateMenu',hMainGui);        
fShow('Image');
fShow('Tracks');
Channel = [Molecule.Channel];
for n = 1:max(Channel)
    k = find(Channel==n);
    for m = 1:length(k)
        X = mean(Molecule(k(m)).Results(:,3)/Molecule(k(m)).PixelSize);
        Y = mean(Molecule(k(m)).Results(:,4)/Molecule(k(m)).PixelSize);
        points{n}(m,:) = double([X Y]);
    end
    if n>1
        idx = 1:size(points{1},1);
        [nidx,D] = knnsearch(points{n},points{1});
        [m,k] = min(D);
        md = median(D);
        nm = 0;
        p = 1;
        ref = [];
        dist = [];
        while ~isempty(m) && (m<md || m<mean(nm)+5*std(nm))
            ref(p,:) = [points{1}(idx(k),1) points{1}(idx(k),2)];
            dist(p,:) = [points{n}(nidx(k),1) points{n}(nidx(k),2)];
            idx(k) = [];
            nidx(k) = [];
            D(k) = [];
            nm(p) = m;
            [m,k] = min(D);
            p = p+1;
        end
        T = estimateGeometricTransform(dist,ref,'similarity');
        hMainGui.Values.TformChannel{n}(:,1:2) = T.T(:,1:2);
        OffsetMap(n-1).Match = [ref dist];
        OffsetMap(n-1).T = T.T;
    end
end
setappdata(hMainGui.fig,'OffsetMap',OffsetMap);
set(hMainGui.Menu.mAlignChannels,'Enable','on');
setappdata(0,'hMainGui',hMainGui);
if strcmp(get(hMainGui.Menu.mShowOffsetMap,'Checked'),'on')
    fShow('OffsetMap',hMainGui);  
end
fShared('UpdateMenu',hMainGui);

function ApplyOffset
global Molecule;
global Filament;
hMainGui = getappdata(0,'hMainGui');
set(hMainGui.Menu.mAlignChannels,'Checked','off');
Molecule = fTransformCoord(Molecule,1,0);
Filament = fTransformCoord(Filament,1,1);
fShared('UpdateMenu',hMainGui);        
fShow('Image');
fShow('Tracks');
for n = 1:length(Molecule)
    c = Molecule(n).Channel;
    Molecule(n).TformMat = hMainGui.Values.TformChannel{c};
end
for n = 1:length(Filament)
    c = Filament(n).Channel;
    Filament(n).TformMat = hMainGui.Values.TformChannel{c};   
end

function Show
hMainGui = getappdata(0,'hMainGui');
if strcmp(get(hMainGui.Menu.mShowOffsetMap,'Checked'),'on')
    set(hMainGui.Menu.mShowOffsetMap,'Checked','Off');
    delete(findobj('Tag','pOffset'));
else
    set(hMainGui.Menu.mShowOffsetMap,'Checked','On');
    fShow('OffsetMap',hMainGui);  
end

function Clear(hMainGui)
OffsetMap = getappdata(hMainGui.fig,'OffsetMap');
Mode = get(gcbo,'UserData');
if strcmp(Mode,'Red')
    OffsetMap.RedXY=[];
elseif strcmp(Mode,'Green')
    OffsetMap.GreenXY=[];
else
    OffsetMap.RedXY=[];    
    OffsetMap.GreenXY=[];    
end
OffsetMap.Match=[];
setappdata(hMainGui.fig,'OffsetMap',OffsetMap);
if strcmp(get(hMainGui.Menu.mShowOffsetMap,'Checked'),'on')
    fShow('OffsetMap',hMainGui);    
end
fShared('UpdateMenu',hMainGui);

function Match(hMainGui)
global Config;
OffsetMap = getappdata(hMainGui.fig,'OffsetMap');
if size(OffsetMap.RedXY,1)~=size(OffsetMap.GreenXY,1)
    fMsgDlg('Number of offset points in the channels do not match','error');
    return;
else
    OffsetMap.Match=[];
    RedXY=OffsetMap.RedXY;
    GreenXY=OffsetMap.GreenXY;
    DiffX=mean(GreenXY(:,1))-mean(RedXY(:,1));
    GreenXY(:,1)=GreenXY(:,1)-DiffX;
    DiffY=mean(GreenXY(:,2))-mean(RedXY(:,2));
    GreenXY(:,2)=GreenXY(:,2)-DiffY;    
    Distance=zeros(size(RedXY,1));
    for i=1:size(RedXY,1)
        Distance(i,:)=sqrt((RedXY(i,1)-GreenXY(:,1)).^2 + (RedXY(i,2)-GreenXY(:,2)).^2);
    end
    while ~isempty(Distance);
        [m,x]=min(min(Distance,[],1));
        [m,y]=min(min(Distance,[],2));   
        OffsetMap.Match=[OffsetMap.Match; RedXY(y,:) GreenXY(x,:)+[DiffX DiffY] ];
        RedXY(y,:)=[]; 
        GreenXY(x,:)=[];
        Distance(y,:)=[];
        Distance(:,x)=[];
    end
    T = estimateGeometricTransform(OffsetMap.Match(:,3:4)/Config.PixSize,OffsetMap.Match(:,1:2)/Config.PixSize,'similarity');
    TformChannel{1} = [1 0 0; 0 1 0; 0 0 1];
    TformChannel{2} = T.T;
    TformChannel{2}(3,3) = 2;
    hMainGui.Values.TformChannel = TformChannel;
    setappdata(0,'hMainGui',hMainGui);
    set(hMainGui.Menu.mAlignChannels,'Enable','on');
end
setappdata(hMainGui.fig,'OffsetMap',OffsetMap);
if strcmp(get(hMainGui.Menu.mShowOffsetMap,'Checked'),'on')
    fShow('OffsetMap',hMainGui);    
end
fShared('UpdateMenu',hMainGui);

function MatchAutomatic(hMainGui)
global Molecule;
OffsetMap = getappdata(hMainGui.fig,'OffsetMap');
if ~isempty(OffsetMap.RedXY)&&~isempty(OffsetMap.GreenXY)
    fMsgDlg('Points present in both channels - clear one','error');
    return;
elseif isempty(OffsetMap.RedXY)&&isempty(OffsetMap.GreenXY)
    fMsgDlg('No points present in any channel - add points to one','error');
    return;
else
    if ~isempty(OffsetMap.RedXY)
        Points = OffsetMap.RedXY;
        mode = 'red';
    else
        Points = OffsetMap.GreenXY;
        mode = 'green';
    end
    for n=1:length(Molecule)
        XY(n,:) = [mean(Molecule(n).Results(:,3)) mean(Molecule(n).Results(:,4))]; 
    end
    for n=1:size(Points,1)
        [dis(n),idx(n)] = min(sqrt( (Points(n,1)-XY(:,1)).^2 + (Points(n,2)-XY(:,2)).^2));
    end
    k = dis > mean(dis)+3*std(dis);
    while any(k)
        Points(k,:)=[];
        dis(k)=[];
        idx(k)=[];
        k = dis > mean(dis)+3*std(dis);
    end
    if strcmp(mode,'red')
        OffsetMap.RedXY = Points;
        OffsetMap.GreenXY = XY(idx,:);
    else
        OffsetMap.RedXY = XY(idx,:);
        OffsetMap.GreenXY = Points;
    end
    OffsetMap.Match=[];
    RedXY=OffsetMap.RedXY;
    GreenXY=OffsetMap.GreenXY;
    DiffX=mean(GreenXY(:,1))-mean(RedXY(:,1));
    GreenXY(:,1)=GreenXY(:,1)-DiffX;
    DiffY=mean(GreenXY(:,2))-mean(RedXY(:,2));
    GreenXY(:,2)=GreenXY(:,2)-DiffY;    
    Distance=zeros(size(RedXY,1));
    for i=1:size(RedXY,1)
        Distance(i,:)=sqrt((RedXY(i,1)-GreenXY(:,1)).^2 + (RedXY(i,2)-GreenXY(:,2)).^2);
    end
    while ~isempty(Distance);
        [m,x]=min(min(Distance,[],1));
        [m,y]=min(min(Distance,[],2));   
        OffsetMap.Match=[OffsetMap.Match; RedXY(y,:) GreenXY(x,:)+[DiffX DiffY] ];
        RedXY(y,:)=[]; 
        GreenXY(x,:)=[];
        Distance(y,:)=[];
        Distance(:,x)=[];
    end
    T = estimateGeometricTransform(ref,dist,'rigid');
end
setappdata(hMainGui.fig,'OffsetMap',OffsetMap);
if strcmp(get(hMainGui.Menu.mShowOffsetMap,'Checked'),'on')
    fShow('OffsetMap',hMainGui);    
end
fShared('UpdateMenu',hMainGui);

function Load(hMainGui)
fRightPanel('CheckOffset',hMainGui);
[FileName, PathName] = uigetfile({'*.mat','FIESTA Offset Map(*.mat)'},'Load FIESTA Offset Map',fShared('GetLoadDir'));
if FileName~=0
    fShared('SetLoadDir',PathName);
    OffsetMap=fLoad([PathName FileName],'OffsetMap');
    setappdata(hMainGui.fig,'OffsetMap',OffsetMap);
    set(hMainGui.Menu.mAlignChannels,'Enable','on');
    if strcmp(get(hMainGui.Menu.mShowOffsetMap,'Checked'),'on')
        fShow('OffsetMap',hMainGui);    
    end
    fShared('UpdateMenu',hMainGui);
end
setappdata(0,'hMainGui',hMainGui);

function Save(hMainGui)
OffsetMap=getappdata(hMainGui.fig,'OffsetMap'); %#ok<NASGU>
[FileName, PathName] = uiputfile({'*.mat','MAT-files (*.mat)'},'Save FIESTA Offset Map',fShared('GetSaveDir'));
if FileName~=0
    fShared('SetSaveDir',PathName);
    file = [PathName FileName];
    if isempty(findstr('.mat',file))
        file = [file '.mat'];
    end
    save(file,'OffsetMap');
end

function newData = CalcNewPos(data,Coeff)
newData(:,1)=data(:,1) * Coeff(1) + data(:,2) * Coeff(2) + Coeff(5);
newData(:,2)=data(:,1) * Coeff(3) + data(:,2) * Coeff(4) + Coeff(6);

function Correct(hMainGui)
global Molecule;
global Filament;
OffsetMap=getappdata(hMainGui.fig,'OffsetMap');
if ~isempty(OffsetMap.Match)
    m=size(OffsetMap.Match,1)*2;
    A=zeros(m,6);
    B=zeros(m,1);
    if strcmp(get(gcbo,'UserData'),'GreenRed')
        A(1:2:m,1:2)=OffsetMap.Match(:,3:4);
        A(2:2:m,3:4)=OffsetMap.Match(:,3:4);    
        B(1:2:m)=OffsetMap.Match(:,1);
        B(2:2:m)=OffsetMap.Match(:,2);
    else
        A(1:2:m,1:2)=OffsetMap.Match(:,1:2);
        A(2:2:m,3:4)=OffsetMap.Match(:,1:2);    
        B(1:2:m)=OffsetMap.Match(:,3);
        B(2:2:m)=OffsetMap.Match(:,4);
    end
    A(1:2:m,5)=1;
    A(2:2:m,6)=1;        
    Coeff=A\B;
    k = find([Molecule.Selected]==1);
    for n = k
        Molecule(n).Results(:,3:4) = CalcNewPos(Molecule(n).Results(:,3:4),Coeff);
    end
    k = find([Filament.Selected]==1);
    for n = k
        Filament(n).Results(:,3:4) = CalcNewPos(Filament(n).Results(:,3:4),Coeff);
        Filament(n).PosStart = CalcNewPos(Filament(n).PosStart,Coeff);
        Filament(n).PosCenter = CalcNewPos(Filament(n).PosCenter,Coeff);
        Filament(n).PosEnd = CalcNewPos(Filament(n).PosEnd,Coeff);
        for m = 1:length(Filament(n).Data)
            Filament(n).Data{m}(:,1:2) =  CalcNewPos(Filament(n).Data{m}(:,1:2),Coeff);
        end
    end
    fShow('Marker',hMainGui,hMainGui.Values.FrameIdx);
    fShow('Tracks');
end