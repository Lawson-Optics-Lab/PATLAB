function WindowLevel(ax,enabled,default,speed_multiplier,legacy_colorbar,enableTagTrigger)
%Enable setting the window level (contrast+brightness) by dragging.
%
% syntax:
% WindowLevel(ax)
% WindowLevel(ax,enabled)
% WindowLevel(ax,enabled,default_caxis_values)
% WindowLevel(ax,enabled,default_caxis_values,speed_multiplier)
% WindowLevel(ax,enabled,default_caxis_values,speed_multiplier,legacy_colorbar)
% WindowLevel(ax,enabled,default_caxis_values,speed_multiplier,legacy_colorbar,enableTagTrigger)
%
% 'ax' and 'enabled' (if provided) are not allowed to be empty.
%
% ax:
%     The axis or array of axes that should have the same caxis. If they are on separate figures,
%     they will be linked, but guidata will fail, so they cannot be used to affect the main window.
%     If ax is an array, ax(1) will be used to load and store guidata in the WindowLevel field.
%     Note: for each figure, only one WindowLevel can be used.
% enabled:
%     A logical or string ('enabled', or 'disabled'). If this is false, any callback will return
%     immediately.
% default_caxis_values:
%     The default window. This should be a 2 element vector containing min_val and max_val. These
%     values control the value for the first and last element of the colormap. If left empty the
%     range of the first element of ax will be used as the default.
% speed_multiplier:
%     A scalar controlling how fast the limits should change. For CT data a value of 10 is advised,
%     the default is 1. If it is not a double, it will be converted to double during input parsing.
% legacy_colorbar:
%     An nx2 cell array. The first col has handles to legacy colorbars, the second col contains the
%     direction ('horiz' or 'vert'). In some older versions (like  ML 6.5), the colorbars don't
%     seem to inherit the color range from their parent axes. Inclusion of these handles will force
%     a range change for those colorbars as well.
% enableTagTrigger:
%     A logical or string ('enabled', or 'disabled'). If this is true, a change in window level
%     will trigger a change in the figure Tag. You can use a listener to catch this event. It will
%     trigger at the release of the mouse button. By default this is false. On releases without
%     addlistener (i.e. pre-R2008a), you can enter a string that will be run with eval. You should
%     only use functions here, but the scope of eval will have access to the handles struct. The
%     string is loaded from the preferences with
%     getpref('HJW','WindowLevel___postset_Tag_evalstr','');.
%
% All persistent data for this function is stored in the WindowLevel field of the guidata of the
% parent figure of the axis (or first axis, in the case of an array of axes).
%
% Because many functions have a tendency to clear the callbacks, it is advisable to re-initialize
% every time after a call to functions like imshow. To prevent this, initialize with imshow and set
% the CData property.
%
% Set [guidata(ax(1))].WindowLevel.enabled to false to disable this function for the entire figure.
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.1.1
% Date:    2020-07-06
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( http://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

%scheme of operation:
%set a flag on ButtonDownFcn and save current cursor position
%use WindowButtonMotionFcn to update the window each time the cursor moves,
%but only if the flag is set
%use WindowButtonUpFcn to reset the flag

if nargin==0 || nargin>6
    error('Incorrect number of inputs. Between 1 and 6 expected.')
end
if nargout~=0
    error('No outputs are supported. Values are saved with guidata, see help.')
end
if numel(ax)==0
    error('The first input (ax) cannot be empty')
end
%Check if all elements of the ax array are handles to axes objects.
for n=1:numel(ax)
    if isprop(ax(n),'Children') && isprop(ax(n),'Type') && ...
            strcmpi(get(ax(n),'Type'),'Axis')
        error('The first input does not contain (a) valid axes handle(s).')
    end
end
if ~exist('enabled','var'),enabled=true;end
if ~islogical(enabled)
    if ischar(enabled)
        enabled=strcmpi(enabled,'enabled')||strcmpi(enabled,'enable');
    else
        error('Unexpected data type for second input (enabled).')
    end
end
if exist('default','var') && ~isempty(default)
    if numel(default)~=2 || ~isnumeric(default)
        error('Third input (default caxis) must be a two element vector.')
    end
else
    default=caxis(ax(1));
end
if exist('speed_multiplier','var') && ~isempty(speed_multiplier)
    if ~isnumeric(speed_multiplier) || numel(speed_multiplier)~=1
        error('Fourth input (speed_multiplier) should be a numeric scalar.')
    end
    speed_multiplier=double(speed_multiplier);
else
    speed_multiplier=1;
end
if exist('legacy_colorbar','var')
    for n=1:size(legacy_colorbar,1)
        if ~isa(legacy_colorbar,'cell')
            error(['Incorrect format for legacy colorbar.'  char(10)...
                'Expected a cell containing directions and handles.' char(10)...
                'See function help text for more information.'])%#ok newline
        end
        if isa(legacy_colorbar{n,2},'double')&&...
                ~strcmpi(get(get(legacy_colorbar{n,2},'Children'),...
                'Tag'),'TMW_COLORBAR')
            error(['Not a legacy colorbar.' char(10)...
                'Use the handle returned by the colorbar function.'])%#ok newline
        end
        if ~ischar(legacy_colorbar{n,1}) || ...
                (strcmpi(legacy_colorbar{n,1},'horiz') &&...
                strcmpi(legacy_colorbar{n,1},'vert'))
            error('Invalid direction for legacy colorbar (use horiz or vert).')
        end
    end
else
    legacy_colorbar=[];
end
if exist('enableTagTrigger','var') && ~isempty(enableTagTrigger)
    if ~islogical(enableTagTrigger)
        if ischar(enableTagTrigger)
            %convert to logical
            enabled=strcmpi(enableTagTrigger,'enabled')...
                || strcmpi(enableTagTrigger,'enable');
        else
            error('Unexpected data type for enableTagTrigger.')
        end
    end
else
    enableTagTrigger=false;
end

%Load or create the handles struct and save persistent settings to it.
handles=guidata(ax(1));
handles.WindowLevel.ax=ax;
handles.WindowLevel.enabled=enabled;
handles.WindowLevel.isBusy=false;%currently processing cursor movement
handles.WindowLevel.listening=false;%listening for cursor movement
handles.WindowLevel.speed_multiplier=speed_multiplier;
handles.WindowLevel.legacy_colorbar=legacy_colorbar;
handles.WindowLevel.enableTagTrigger=enableTagTrigger;

handles.WindowLevel.use_legacy=ifversion('<',7,'Octave','<=',3);%use legacy colorbar
%addlistener was introduced in R2008a, so before that eval a string
handles.WindowLevel.UseEval_TagTrigger=ifversion('<','R2008a','Octave','<=',3);

%Set the button down function for all axes in the list and their children.
%No other objects in the figure should activate the windowing.
for m=1:numel(ax)
    caxis(ax(m),default);
    set(ax(m),'ButtonDownFcn',@WindowLevelButtonDown)
    children_list=get(ax(m),'Children');
    for n=1:length(children_list)
        set(children_list(n),'ButtonDownFcn',@WindowLevelButtonDown)
    end
end

%get the parent fig of ax(1) and set the callbacks
parentfig=ax(1);
while ~strcmpi(get(parentfig,'Type'),'Figure')
    parentfig=get(parentfig,'Parent');
end
set(parentfig,'WindowButtonUpFcn',@WindowLevelButtonUp)
set(parentfig,'WindowButtonMotionFcn',@WindowLevelButtonMotion)
guidata(ax(1),handles)
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
% tf=ifversion(test,Rxxxxab)
% tf=ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Output:
% tf       - If the current version satisfies the test this returns true.
%            This works similar to verLessThan.
%
% Inputs:
% Rxxxxab - Char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the
%           numeric version.
% test    - Char array containing a logical test. The interpretation of this is equivalent to
%           eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',9.8) returns true only when run on R2020a
% ifversion('<',0,'Octave','>',0) returns true only on Octave
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.0.2
% Date:    2020-05-20
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

%The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
%This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
%remove the potential for float rounding errors.
%Store in persistent for fast recall (don't use getpref, as that is slower than generating the
%variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    %test if Octave is used instead of Matlab
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    %get current version number
    v_num=version;
    ii=strfind(v_num,'.');
    if numel(ii)~=1,v_num(ii(2):end)='';ii=ii(1);end
    v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
    v_num=v_num(1)+v_num(2)/100;
    v_num=round(100*v_num);%remove float rounding errors
    
    %get dictionary to use for ismember
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;'R14SP3' 701;...
        'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;'R2008a' 706;'R2008b' 707;...
        'R2009a' 708;'R2009b' 709;'R2010a' 710;'R2010b' 711;'R2011a' 712;'R2011b' 713;...
        'R2012a' 714;'R2012b' 800;'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;...
        'R2015a' 805;'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;...
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        %Use the same test as for Matlab, which will probably fail.
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    elseif nargin==4
        %undocumented shorthand syntax: skip the 'Octave' argument
        [test,v]=deal(Oct_flag,Oct_test);
        %convert 4.1 to 401
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        %convert 4.1 to 401
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    %convert R notation to numeric and convert 9.1 to 901
    if isnumeric(Rxxxxab)
        v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
    else
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    end
end
switch test
    case '=='
        tf= v_num == v;
    case '<'
        tf= v_num <  v;
    case '<='
        tf= v_num <= v;
    case '>'
        tf= v_num >  v;
    case '>='
        tf= v_num >= v;
end
end
function WindowLevelButtonDown(src,unused_evendata)%#ok<INUSD> ~ ML6.5
handles=guidata(src);
%Check isfield in case the guidata is wiped, but the callback isn't.
if ~isfield(handles,'WindowLevel'),return,end
enabled=handles.WindowLevel.enabled;
if ~islogical(enabled)
    enabled=strcmpi(enabled,'enabled')||strcmpi(enabled,'enable');
end
if ~enabled
    return
end

%start listening for cursor motion
handles.WindowLevel.listening=true;

%set units to pixels, but save the old units
handles.WindowLevel.oldFigureUnits=get(gcf,'Units');
set(gcf,'Units','Pixels');

%get the current cursor position and current display range
handles.WindowLevel.currentStartPointOfMotion_pos=get(gcf,'CurrentPoint');
handles.WindowLevel.currentStartPointOfMotion_caxis=caxis(gca);
guidata(src,handles)
end
function WindowLevelButtonMotion(src,unused_evendata)%#ok<INUSD> ~ ML6.5
handles=guidata(src);
%Check isfield in case the guidata is wiped, but the callback isn't.
enabled=isfield(handles,'WindowLevel') && handles.WindowLevel.enabled;
if ~islogical(enabled)
    enabled=strcmpi(enabled,'enabled')||strcmpi(enabled,'enable');
end
if      ~enabled ||... motion detection not activated
        handles.WindowLevel.isBusy ||...%callback currently being executed
        ~handles.WindowLevel.listening %listening for motion
    return
end
handles.WindowLevel.isBusy=true;
guidata(src,handles)%prevent too many callbacks to be executed
% try
newpos=get(gcf,'CurrentPoint');
oldpos=handles.WindowLevel.currentStartPointOfMotion_pos;
old_lim=handles.WindowLevel.currentStartPointOfMotion_caxis;
deltapos=double(handles.WindowLevel.speed_multiplier)*(oldpos-newpos);
level=mean(old_lim);
window=diff(old_lim);

% event |  intended result   | delta position
%---------------------------------------------
% down  |  lower the level   |  deltapos(2)>0
% left  |  lower the window  |  deltapos(1)>0

level=level-deltapos(2);
window=window-deltapos(1);
if window<0,window=abs(window);end%make sure window is positive
if window==0,window=2*eps*level;end%make sure new_lim is a valid vector
new_lim=level+[-0.5 0.5]*window;
for n=1:numel(handles.WindowLevel.ax)
    caxis(handles.WindowLevel.ax(n),new_lim);
end

%Check if the version is old enough to make legacy colorbar handling plausibly necessary. It is
%difficult to find when the behavior changed substantially, but it is somewhere between ML6.5 and
%R2011a, my guesstimate is ML7. The method below is unreliable and needs fixing for colorbars that
%*do* work correctly, so here, the entire section will be ignored for R14 (v7.0) and later. (don't
%use this code for Octave either, as that doesn't seem to suffer from this (only Octave 4.2.1 was
%tested))
lc=handles.WindowLevel.legacy_colorbar;
if ~isempty(lc) && handles.WindowLevel.use_legacy
    new_lim__=linspace(new_lim(1),new_lim(2),4);
    new_lim__(1)=ceil(new_lim__(1));new_lim__=floor(new_lim__);
    new_lim__=unique(new_lim__);
    for n=1:size(lc,1)
        if strcmpi(lc{n,1},'horiz'),d='X';else,d='Y';end
        val_range=get(lc{n,2},[d 'Lim']);
        new_lim_=(new_lim__-new_lim(1))/window*diff(val_range)+min(val_range);
        set(lc{n,2},[d 'Tick'],new_lim_)
        set(lc{n,2},[d 'TickLabel'],num2str(new_lim__(:)))
    end
end
% catch
%     warning('HJW:WindowLevel:UnkownErrorWarning',...
%         'An error occurred, but no error handling is implemented yet.')
% end
handles.WindowLevel.isBusy=false;%release again for a new motion callback
guidata(src,handles)
end
function WindowLevelButtonUp(src,unused_evendata)%#ok<INUSD> ~ ML6.5
handles=guidata(src);
%Check isfield in case the guidata is wiped, but the callback isn't.
if ~isfield(handles,'WindowLevel'),return,end
enabled=handles.WindowLevel.enabled;
if ~islogical(enabled)
    enabled=strcmpi(enabled,'enabled')||strcmpi(enabled,'enable');
end
if ~enabled
    return
end

%stop listening for cursor motion
handles.WindowLevel.listening=false;
guidata(src,handles)

try
    %reset the figure Units
    set(gcf,'Units',handles.DragNav.oldFigureUnits);
catch
end
if handles.WindowLevel.enableTagTrigger
    %trigger a change in the tag for when the listener is active addlistener was introduced in
    %R2008a, so before that eval a string
    if handles.WindowLevel.UseEval_TagTrigger
        evalstr=getpref('HJW',...%toolbox name or author initials as group ID
            ['WindowLevel___',...%function group
            'evalstr'],...%preference name
            []);
        if ~isempty(evalstr)
            eval(evalstr);
        end
    else
        %randchar=randi([65 122],1,40);%randi was introduced in R2008b
        randchar=uint8(rand(1,40)*(122-65)+65);
        randchar=char(randchar(randchar<91|randchar>96));%select only valid letters
        set(gcf,'Tag',randchar)
    end
end
end