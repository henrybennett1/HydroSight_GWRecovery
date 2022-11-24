%function [f, status] = HydroSight()
function varargout = HydroSight(doingTesting)

    if nargin==0
        doingTesting=false;
    end
    % Check Matlab is 2018a or later. 
    if ~isdeployed
        v=version();
        ind = strfind(v,'.');
        v_major = str2double(v(1:(ind(1)-1)));
        v_minor = str2double(v((ind(1)+1):(ind(2)-1)));
        if v_major<9 && v_minor<4 %ie 2018a;
            errordlg('HydroSight requires Matlab 2018a or later.','Please update Matlab');
            return;
        end
    end

    % Add Paths
    addpath(genpath([pwd, filesep, 'algorithms']));
    addpath(genpath([pwd, filesep, 'GUI']));
        
    % Load GUI
    try
        % Use GUI Layout Toolbox if it exists.
        if ~isdeployed && ~isempty(ver('layout'))           
            rmpath(genpath(fullfile( pwd, 'GUI','GUI Layout Toolbox 2.3.4')));
        end
        f = HydroSight_GUI();
        if doingTesting
            varargout = {f};
        end

    catch ME
        % Check the toolbox for GUIs exists
        if ~isdeployed && isempty(ver('layout')) 
            msgbox({'The GUI cannot be opened. This may be because the GUI Layout Toolbox is ', ...
                    'not installed within Matlab. Please download it using the URL below and', ...
                    'then install it and re-start HydroSight.', ...
                    '', ...
                    'https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox'}, ...
                    'Toolbox missing: gui-layout-toolbox', 'error');

        else
            if size(ME.stack,1)>1
                functionName = ME.stack(end-1,1).name;
                functionLine = num2str(ME.stack(end-1,1).line);
            else
                functionName='';
                functionLine ='';
            end
            msgbox({'An unexpected error occured. To help with fixing the issue please', ...
                    'copy the error message below and submit a bug report to:', ...
                    'https://github.com/peterson-tim-j/HydroSight/issues', ...
                    '', ...
                    ['Message:', ME.message], ...
                    '', ...
                    ['Function:', functionName], ...
                    '', ...
                    ['Line Number:', functionLine]}, ...
                    'Unknown error', 'error');
            
        end
        if doingTesting
            varargout = {};
        end
    end

    % Check if the required toolboxes are installed.
    if ~isdeployed
        v = ver;
        [installedToolboxes{1:length(v)}] = deal(v.Name);

        if ~any(strcmp(installedToolboxes, 'Statistics and Machine Learning Toolbox'))
            warndlg({'Statistics and Machine Learning Toolbox is required.', ...
                'This toolbox is required for the most model algorithms.'},'Statistics Toolbox not installed','modal')
        end

        if ~any(strcmp(installedToolboxes, 'Parallel Computing Toolbox'))
            warndlg({'Parallel Computing Toolbox is suggested.', ...
                'This toolbox is used to reduce the calibration time for the most models.'},'Parallel Toolbox not installed','modal')
        end

        if ~any(strcmp(installedToolboxes, 'Curve Fitting Toolbox'))
            warndlg({'Curve Fitting Toolbox is required.', ...
                'This toolbox is required for the outlier detection algorithm.'},'Curve Fitting Toolbox not installed','modal')
        end

    end

end
