function results = testHydroSight(runTests)
% Setup test environment and run suite.
%
% Setup includes closig the parpool and ensuring it does not autostart. The
% GUI layout toolbox is also installed. If the setup is successful, and
% only the setup is being done, then result equals 0, else -1;
%
% If runTests==true, then the test suite is also built and run. If so, then
% result equals a cell array of the TestResult objects.
%
    disp('SETTING UP TEST ENVIONMENT');
    disp('--------------------------')
    
    results = 0;
    
    try
        disp('Import testing plugins ...')
        import matlab.unittest.plugins.CodeCoveragePlugin %#ok<SIMPT> 
    
        % Close parallel pool. This is done to ensure code within a parfor is
        % included in code coverage report
        disp('Closing parpool...')
        poolobj = gcp('nocreate');
        delete(poolobj);
    
        % Ensure parpool does not auto-start
        disp('Ensure parpool is not automatically started ...')
        ps = parallel.Settings;
        ps_Auto = ps.Pool.AutoCreate;
        ps.Pool.AutoCreate = false;
    
        % Add paths to unit tests
        disp('Adding testing to path ...')
        addpath('testing');
        currentPath = pwd();
    
        % Install GUI layout if not installed
        if isempty(ver('layout'))
            disp('Installng GUI Layout toolbox ...')
            matlab.addons.toolbox.installToolbox(fullfile(currentPath,'testing','GUI Layout Toolbox 2.3.5.mltbx'),true);
        end
        disp('SETUP COMPLETE.')
    catch ME
        results = -1;
        disp(['Test setup failed with error: ',ME.message]);
    end
    disp('--------------------------')
    disp(' ')

    if runTests || nargin==0
        disp('STARTING TESTS');
        disp('--------------------------')
               
        % Setyp tests
        disp('Making test suite ...')
        suite = testsuite({'runOutlierDetection','buildExampleModels','simulateExampleModels','calibrateNewModel'});
        %suite = testsuite({'calibrateNewModel'});        
        runner = testrunner("textoutput");

        % add coverage report
        disp('Adding coverage report ...')
        runner.addPlugin(CodeCoveragePlugin.forFolder({'GUI','algorithms'},'IncludeSubfolders',true));

        disp('Run tests...')
        results = runner.run(suite);
    
        disp('FINSIHED TESTS');
        disp('--------------------------')
        disp('')
        
        % Rest file paths.
        disp('Restoring paths ...')
        rmpath(fullfile(currentPath,'testing'));
        cd(currentPath);
    
        % Restore parpool settings
        disp('Restoring parpool settings ...')
        ps.Pool.AutoCreate = ps_Auto;
    end
end