classdef Intestinal_Calcium_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        IntestinalCalciumAppUIFigure   matlab.ui.Figure
        FileMenu                       matlab.ui.container.Menu
        saveSettings                   matlab.ui.container.Menu
        savedefaultsettingsMenu        matlab.ui.container.Menu
        loadSettings                   matlab.ui.container.Menu
        TabGroup                       matlab.ui.container.TabGroup
        TrackingTab                    matlab.ui.container.Tab
        AnalyzeOASdataCheckBox         matlab.ui.control.CheckBox
        StartFrameSpinner              matlab.ui.control.Spinner
        StartFrameSpinnerLabel         matlab.ui.control.Label
        RunAnalysisButton              matlab.ui.control.Button
        StartFileSpinner               matlab.ui.control.Spinner
        StartFileSpinnerLabel          matlab.ui.control.Label
        FlatfieldCorrectionEditFieldLabel  matlab.ui.control.Label
        FlatfieldCorrectionEditField   matlab.ui.control.NumericEditField
        InputFramerateEditFieldLabel   matlab.ui.control.Label
        InputFramerateEditField        matlab.ui.control.NumericEditField
        CropPixelsEditFieldLabel       matlab.ui.control.Label
        CropPixelsEditField            matlab.ui.control.NumericEditField
        ThresholdingButtonGroup        matlab.ui.container.ButtonGroup
        AdaptiveThresholdButton        matlab.ui.control.RadioButton
        AutoThresholdButton            matlab.ui.control.RadioButton
        DisplayOptionsPanel            matlab.ui.container.Panel
        VideoFramerateEditFieldLabel   matlab.ui.control.Label
        VideoFramerateEditField        matlab.ui.control.NumericEditField
        RecordVideoCheckBox            matlab.ui.control.CheckBox
        TroubleshootCheckBox           matlab.ui.control.CheckBox
        ShowNormalsCheckBox            matlab.ui.control.CheckBox
        ShowAxialSignalCheckBox        matlab.ui.control.CheckBox
        PlotResultsCheckBox            matlab.ui.control.CheckBox
        UploadResultsCheckBox          matlab.ui.control.CheckBox
        remoteDir                      matlab.ui.control.EditField
        setRemoteDirButton             matlab.ui.control.Button
        LoadTiffCheckBox               matlab.ui.control.CheckBox
        IsRemoteCheckBox               matlab.ui.control.CheckBox
        tiffDir                        matlab.ui.control.EditField
        setTiffDirButton               matlab.ui.control.Button
        PlottingTab                    matlab.ui.container.Tab
        Panel_2                        matlab.ui.container.Panel
        PrefixEditField                matlab.ui.control.EditField
        PrefixEditFieldLabel           matlab.ui.control.Label
        numColumns                     matlab.ui.control.NumericEditField
        ColumnsEditFieldLabel          matlab.ui.control.Label
        mutantgenotypesEditField       matlab.ui.control.EditField
        mutantgenotypesEditFieldLabel  matlab.ui.control.Label
        PlotMultipleGenotypesButton    matlab.ui.control.Button
        Panel                          matlab.ui.container.Panel
        windowsizeEditField            matlab.ui.control.NumericEditField
        windowsizeLabel                matlab.ui.control.Label
        spikeEditField                 matlab.ui.control.NumericEditField
        spikeEditFieldLabel            matlab.ui.control.Label
        PlotSingleTraceButton          matlab.ui.control.Button
        CombinewormDataButton          matlab.ui.control.Button
        MergeControlButton             matlab.ui.control.Button
        PlotMatchedControlButton       matlab.ui.control.Button
        controlnameEditField           matlab.ui.control.EditField
        controlnameEditFieldLabel      matlab.ui.control.Label
        setOutputDir                   matlab.ui.control.Button
        outputDir                      matlab.ui.control.EditField
        DataDir                        matlab.ui.control.EditField
        setDataDir                     matlab.ui.control.Button
        PlotSettingsTab                matlab.ui.container.Tab
        AxialSignalColorDropDown       matlab.ui.control.DropDown
        AxialSignalColorLabel          matlab.ui.control.Label
        EqualizeExpDurationCheckBox    matlab.ui.control.CheckBox
        AnalyzeOASdataCheckBox_2       matlab.ui.control.CheckBox
        PlotMultipleGenotypesButton_2  matlab.ui.control.Button
        PlotMatchedControlButton_2     matlab.ui.control.Button
        numPlotsEditField              matlab.ui.control.NumericEditField
        PlotsEditFieldLabel            matlab.ui.control.Label
        BulkSignalPlotsButtonGroup     matlab.ui.container.ButtonGroup
        OverlayButton                  matlab.ui.control.RadioButton
        SeparateButton                 matlab.ui.control.RadioButton
        SortDirectionButtonGroup       matlab.ui.container.ButtonGroup
        shuffleButton                  matlab.ui.control.RadioButton
        ascendButton                   matlab.ui.control.RadioButton
        descendButton                  matlab.ui.control.RadioButton
        FractiontoQuerrylengthLabel    matlab.ui.control.Label
        toQuerry                       matlab.ui.control.NumericEditField
        AutoFixAxialSignalCheckBox     matlab.ui.control.CheckBox
        NormalizeCheckBox              matlab.ui.control.CheckBox
        FrameRateEditField             matlab.ui.control.NumericEditField
        FrameRateEditFieldLabel        matlab.ui.control.Label
        PeakDetectionPanel             matlab.ui.container.Panel
        spikeProfileWindow             matlab.ui.control.NumericEditField
        PlotwindowsecLabel             matlab.ui.control.Label
        peakwidth                      matlab.ui.control.NumericEditField
        MinWidthframesEditFieldLabel   matlab.ui.control.Label
        peakdistance                   matlab.ui.control.NumericEditField
        MinDistanceframesEditFieldLabel  matlab.ui.control.Label
        peakthreshold                  matlab.ui.control.NumericEditField
        PeakProminenceAULabel          matlab.ui.control.Label
        HistogramBinsPanel             matlab.ui.container.Panel
        IncrementEditField             matlab.ui.control.NumericEditField
        IncrementEditFieldLabel        matlab.ui.control.Label
        MaxEditField                   matlab.ui.control.NumericEditField
        MaxEditFieldLabel              matlab.ui.control.Label
        MinEditField                   matlab.ui.control.NumericEditField
        MinEditFieldLabel              matlab.ui.control.Label
        SortTracesByButtonGroup        matlab.ui.container.ButtonGroup
        sortAmp                        matlab.ui.control.RadioButton
        sortFreq                       matlab.ui.control.RadioButton
        dontSort                       matlab.ui.control.RadioButton
        XAxesPanel                     matlab.ui.container.Panel
        axialXTickInt                  matlab.ui.control.EditField
        BulkSignallimitsLabel_2        matlab.ui.control.Label
        bulkXLim                       matlab.ui.control.EditField
        BulkSignallimitsLabel          matlab.ui.control.Label
        YAxesPanel                     matlab.ui.container.Panel
        bulkYLim                       matlab.ui.control.EditField
        BulkSignallimsLabel            matlab.ui.control.Label
        axialYLim                      matlab.ui.control.EditField
        AxialSignallimLabel            matlab.ui.control.Label
        MultigenotypefigurepositionPanel  matlab.ui.container.Panel
        setTracePos                    matlab.ui.control.Button
        setGraphPos                    matlab.ui.control.Button
        TracesEditField                matlab.ui.control.EditField
        TracesEditFieldLabel           matlab.ui.control.Label
        GraphsEditField                matlab.ui.control.EditField
        GraphsEditFieldLabel           matlab.ui.control.Label
        miscsettingsTab                matlab.ui.container.Tab
        SpikeKineticsPanel             matlab.ui.container.Panel
        validateRiseFall               matlab.ui.control.CheckBox
        validatePropagationRate        matlab.ui.control.CheckBox
        showFitParams                  matlab.ui.control.CheckBox
        validateRiseFallButton         matlab.ui.control.Button
        ColorsPanel                    matlab.ui.container.Panel
        wtColor                        matlab.ui.control.EditField
        wildtypeedgecolorEditFieldLabel  matlab.ui.control.Label
        wtEdgeColor                    matlab.ui.control.EditField
        wildtypecolorEditFieldLabel    matlab.ui.control.Label
        mtColor                        matlab.ui.control.EditField
        mutantcolorEditFieldLabel      matlab.ui.control.Label
    end

    
    properties (Access = private)
        defaultDataDir='Z:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants'; % Description
        defaultOutputDir = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants';
    end

    methods (Access = public)

        function parsedInputs = parseInputs(app)

            parsedInputs.tiffDir = app.tiffDir.Value;

            parsedInputs.uploadResults = app.UploadResultsCheckBox.Value;

            if ~isfolder(app.remoteDir.Value)
                parsedInputs.uploadResults = 0;
                parsedInputs.remoteDir = [];
            else
                parsedInputs.remoteDir = app.remoteDir.Value;
            end

            parsedInputs.loadTiff = app.LoadTiffCheckBox.Value;
            parsedInputs.isRemote = app.IsRemoteCheckBox.Value;
            parsedInputs.crop = app.CropPixelsEditField.Value;
            parsedInputs.inputFramerate = app.InputFramerateEditField.Value;
            parsedInputs.autoThreshold = app.AutoThresholdButton.Value;
            parsedInputs.adaptiveThreshold = app.AdaptiveThresholdButton.Value;
            parsedInputs.flatField = app.FlatfieldCorrectionEditField.Value;
            parsedInputs.plotResults = app.PlotResultsCheckBox.Value;
            parsedInputs.showAxialSignal = app.ShowAxialSignalCheckBox.Value;
            parsedInputs.showNormals = app.ShowNormalsCheckBox.Value;
            parsedInputs.troubleshoot = app.TroubleshootCheckBox.Value;
            parsedInputs.recordVideo = app.RecordVideoCheckBox.Value;
            parsedInputs.videoFramerate = app.VideoFramerateEditField.Value;
            parsedInputs.startFrame = app.StartFrameSpinner.Value;
            parsedInputs.startFile = app.StartFileSpinner.Value;
            parsedInputs.isOAS = app.AnalyzeOASdataCheckBox.Value;




        end
    end

    methods (Access = private)

        function plotSettings = parsePlotSettings(app)

            plotSettings.axylimit = str2num(app.axialYLim.Value);
            plotSettings.traceylimit = str2num(app.bulkYLim.Value);


            plotSettings.xlimits = str2num(app.bulkXLim.Value);
            plotSettings.axialXticint = str2double(app.axialXTickInt.Value);

            plotSettings.peakthreshold = app.peakthreshold.Value;
            plotSettings.peakdistance = app.peakdistance.Value;
            plotSettings.peakwidth = app.peakwidth.Value;
            plotSettings.spikeProfileWindow = app.spikeProfileWindow.Value;

            plotSettings.normalize = app.NormalizeCheckBox.Value;

            plotSettings.autoFixAxialSignal = app.AutoFixAxialSignalCheckBox.Value;
            plotSettings.axSigToQuerry = app.toQuerry.Value;

            plotSettings.framerate = app.FrameRateEditField.Value;

            %% plot color 
            plotSettings.wtcolor = str2num(app.wtColor.Value);
            plotSettings.mtcolor = str2num(app.mtColor.Value);
            plotSettings.mtedgecolor = str2num(app.wtEdgeColor.Value);

            plotSettings.axSigCMap = app.AxialSignalColorDropDown.Value;
            


            plotSettings.binedges = app.MinEditField.Value:app.IncrementEditField.Value:app.MaxEditField.Value;
            %% sorting

            if app.dontSort.Value
                plotSettings.sortType = 0;
            elseif app.sortFreq.Value
                plotSettings.sortType = 1;
            elseif app.sortAmp.Value
                plotSettings.sortType = 2;
            end

            if app.descendButton.Value
                plotSettings.sortDir = 'descend';
            elseif app.ascendButton.Value
                plotSettings.sortDir = 'ascend';
            elseif app.shuffleButton.Value
                plotSettings.sortDir = 'shuffle';
            end

            if app.SeparateButton.Value
                plotSettings.overlayplots = 0;
            elseif app.OverlayButton.Value
                plotSettings.overlayplots = 1;
            end

            plotSettings.tolimit = app.numPlotsEditField.Value;
            plotSettings.controlname = app.controlnameEditField.Value;
            
            %% Genotype Parsing for plot_MultiGenotype
            gtyps = split(app.mutantgenotypesEditField.Value, {' ' , ','});
            gtyps = gtyps(~cellfun(@isempty,gtyps));

            if length(gtyps) == 1
                gtyps = gtyps{:};
            end

            plotSettings.genotypes = gtyps;
            plotSettings.prefix = app.PrefixEditField.Value;

            %% Layout Settings for plot_MultiGenotype
            if app.numColumns.Value == 0 
                numcolumns = length(gtyps);
            else
                numcolumns = app.numColumns.Value;
            end
            plotSettings.numColumns = numcolumns;

            plotSettings.graphPos = str2num(app.GraphsEditField.Value);
            plotSettings.tracePos = str2num(app.TracesEditField.Value);

            %% OAS related settings
            plotSettings.isOAS = app.AnalyzeOASdataCheckBox.Value;
            plotSettings.trimExperimentLength = app.EqualizeExpDurationCheckBox.Value;

            %% Spike Kinetics and Validation Settings
            plotSettings.validatePropagationRate = app.validatePropagationRate.Value;
            plotSettings.validateRiseFall = app.validateRiseFall.Value;
            plotSettings.showFitParams = app.showFitParams.Value;

            %% single trace plotting
            plotSettings.singleSpike = app.spikeEditField.Value;
            plotSettings.spikeWindow = app.windowsizeEditField.Value;


        end

            
        function applyPlotSettings(app,plotSettings)
            %% axis limits
            app.axialYLim.Value=  num2str(plotSettings.axylimit);
            app.bulkYLim.Value= num2str(plotSettings.traceylimit);
            app.bulkXLim.Value = num2str(plotSettings.xlimits);
            app.axialXTickInt.Value = num2str(plotSettings.axialXticint);
            %% peak detection
            app.peakthreshold.Value = plotSettings.peakthreshold;
            app.peakdistance.Value= plotSettings.peakdistance;
            app.peakwidth.Value= plotSettings.peakwidth;
            app.spikeProfileWindow.Value= plotSettings.spikeProfileWindow;

            app.NormalizeCheckBox.Value = plotSettings.normalize;

            app.AutoFixAxialSignalCheckBox.Value= plotSettings.autoFixAxialSignal;
            app.toQuerry.Value = plotSettings.axSigToQuerry;

            app.FrameRateEditField.Value = plotSettings.framerate;
            %% plot color
            app.wtColor.Value = num2str(plotSettings.wtcolor);
            app.mtColor.Value= num2str(plotSettings.mtcolor);
            app.wtEdgeColor.Value = num2str(plotSettings.mtedgecolor);

            app.AxialSignalColorDropDown.Value = plotSettings.axSigCMap;
            
            %% histogram Bins
            app.MinEditField.Value = plotSettings.binedges(1);
            app.MaxEditField.Value = plotSettings.binedges(end);
            app.IncrementEditField.Value = plotSettings.binedges(end)/(length(plotSettings.binedges)-1);
            
            %% sorting

            switch plotSettings.sortType
                case 0
                    app.dontSort.Value = 1;
                case 1
                    app.sortFreq.Value = 1;
                case 2
                    app.sortAmp.Value = 1;
            end


            if strcmpi(plotSettings.sortDir,'descend')
                app.descendButton.Value = 1;
            elseif strcmpi(plotSettings.sortDir,'ascend')
                app.ascendButton.Value = 1;
            elseif strcmpi(plotSettings.sortDir,'shuffle')
                app.shuffleButton.Value = 1;
            end
            
            if plotSettings.overlayplots
                app.OverlayButton.Value = 1;
            else
                app.SeparateButton.Value = 1;
            end

            app.numPlotsEditField.Value = plotSettings.tolimit;
            app.controlnameEditField.Value = plotSettings.controlname;


            app.GraphsEditField.Value = num2str(plotSettings.graphPos);
            app.TracesEditField.Value = num2str(plotSettings.tracePos);

            %% OAS related settings
            app.AnalyzeOASdataCheckBox.Value = plotSettings.isOAS;
            app.EqualizeExpDurationCheckBox.Value = plotSettings.trimExperimentLength;
            

            
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            [path, ~, ~] = fileparts(mfilename("fullpath"));
            load(fullfile(path, 'defaultPlotSettings.mat'), 'defaultPlotSettings');
            applyPlotSettings(app, defaultPlotSettings)
        end

        % Button pushed function: setTiffDirButton
        function setTiffDirButtonButtonPushed(app, event)
            prevDir = app.tiffDir.Value;
            if isfolder(prevDir)
                searchStart = prevDir;
            else
                searchStart = app.defaultOutputDir; % 'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants';
            end

            selectedDir = uigetdir(searchStart, 'Select tiff directory');
            figure(app.IntestinalCalciumAppUIFigure)

            if isfolder(selectedDir)
                app.tiffDir.Value = selectedDir;
            end


        end

        % Button pushed function: setRemoteDirButton
        function setRemoteDirButtonButtonPushed(app, event)
            prevDir = app.remoteDir.Value;
            if isfolder(prevDir)
                searchStart = prevDir;
            else
                searchStart = 'Z:\';
            end

            selectedDir = uigetdir(searchStart, 'Select tiff directory');
            figure(app.IntestinalCalciumAppUIFigure)

            if selectedDir ~= 0
                app.remoteDir.Value = selectedDir;
            end

        end

        % Button pushed function: RunAnalysisButton
        function RunAnalysisButtonPushed(app, event)
            parsedInputs = parseInputs(app);
            if parsedInputs.isOAS == 0
            freelyMovingAnalysis_Func(parsedInputs)
            elseif parsedInputs.isOAS==1
                OAS_Analysis_Func(parsedInputs)
            end


        end

        % Button pushed function: setDataDir
        function setDataDirButtonPushed(app, event)
            prevDir = app.DataDir.Value;
            if isfolder(prevDir)
                searchStart = prevDir;
            else
                searchStart = app.defaultDataDir; %'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants';
            end

            selectedDir = uigetdir(searchStart, 'Select wormdata.mat directory');
            figure(app.IntestinalCalciumAppUIFigure)

            if isfolder(selectedDir)
                app.DataDir.Value = selectedDir;
            end
        end

        % Button pushed function: setOutputDir
        function setOutputDirButtonPushed(app, event)
            prevDir = app.outputDir.Value;
            if isfolder(prevDir)
                searchStart = prevDir;
            else
                searchStart = app.defaultOutputDir;% 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\';
            end

            selectedDir = uigetdir(searchStart, 'Select mergedData directory');
            figure(app.IntestinalCalciumAppUIFigure)

            if isfolder(selectedDir)
                app.outputDir.Value = selectedDir;
            end

        end

        % Button pushed function: CombinewormDataButton
        function CombinewormDataButtonPushed(app, event)
            datadir = app.DataDir.Value;
            if ~isfolder(datadir)
                datadir = uigetdir(app.defaultDataDir,...
                    'Oops! looks like you forgot to select your data directory!');
                figure(app.IntestinalCalciumAppUIFigure)
            end
            

            outputdir = app.outputDir.Value;

            if ~isfolder(outputdir)
                outputdir = [];
            end


            combine_Wormdata(datadir,outputdir,app.controlnameEditField.Value);

            app.outputDir.Value = outputdir;
            app.DataDir.Value = datadir;
        end

        % Button pushed function: MergeControlButton
        function MergeControlButtonPushed(app, event)
            outputdir = app.outputDir.Value;
            if ~isfolder(outputdir)
                if isfolder(app.DataDir.Value)
                    outputdir = app.DataDir.Value;
                else
                    outputdir = uigetdir(app.defaultDataDir,...
                        'Oops! looks like you forgot to select your data directory!');
                    app.DataDir.Value = outputdir;
                    figure(app.IntestinalCalciumAppUIFigure)
                end
            end
            mergeControl(outputdir, app.controlnameEditField.Value)
        end

        % Value changed function: NormalizeCheckBox
        function NormalizeCheckBoxValueChanged(app, event)
            value = app.NormalizeCheckBox.Value;
            if value == 0
                app.bulkYLim.Value = '-500 5000';
                app.peakthreshold.Value = 500;
            elseif value == 1
                app.bulkYLim.Value = '-0.2 1';
                app.peakthreshold.Value = 0.1;
            end



        end

        % Button pushed function: PlotMatchedControlButton, 
        % ...and 1 other component
        function PlotMatchedControlButtonPushed(app, event)
            outputdir = app.outputDir.Value;
            if ~isfolder(outputdir)
                if isfolder(app.DataDir.Value)
                    outputdir = app.DataDir.Value;
                else
                    outputdir = uigetdir(app.defaultDataDir,...
                        'Oops! looks like you forgot to select your data directory!');
                    app.DataDir.Value = outputdir;
                    figure(app.IntestinalCalciumAppUIFigure)
                end
            end
            settings = parsePlotSettings(app);
            plot_MatchedControl(app.outputDir.Value,settings);

        end

        % Button pushed function: PlotMultipleGenotypesButton, 
        % ...and 1 other component
        function PlotMultipleGenotypesButtonPushed(app, event)
            
            outputdir = app.outputDir.Value;
            if ~isfolder(outputdir)
                if isfolder(app.DataDir.Value)
                    outputdir = app.DataDir.Value;
                else
                    outputdir = uigetdir(app.defaultDataDir,...
                        'Oops! looks like you forgot to select your data directory!');
                    app.DataDir.Value = outputdir;
                    figure(app.IntestinalCalciumAppUIFigure)
                end
            end
            settings = parsePlotSettings(app);
            plot_MultiGenotype(outputdir, settings)


        end

        % Value changed function: AnalyzeOASdataCheckBox, 
        % ...and 1 other component
        function AnalyzeOASdataCheckBox_2ValueChanged(app, event)
            OASvalue = app.AnalyzeOASdataCheckBox_2.Value;
            normalizeValue = app.NormalizeCheckBox.Value;
            if OASvalue == 1
                app.bulkYLim.Value = '0 15';
                app.axialYLim.Value = '0 45';
                app.peakthreshold.Value = 2;
                app.EqualizeExpDurationCheckBox.Value = 1;
            end

            if OASvalue == 0
                app.EqualizeExpDurationCheckBox.Value = 0;
                if normalizeValue == 0
                    app.bulkYLim.Value = '-500 5000';
                    app.axialYLim.Value = '-500 12000';
                    app.peakthreshold.Value = 500;
                    
                elseif normalizeValue == 1
                    app.bulkYLim.Value = '-0.2 1';
                    app.peakthreshold.Value = 0.1;
                end
            end
            
        end

        % Button pushed function: validateRiseFallButton
        function validateRiseFallButtonPushed(app, event)
            app.validateRiseFall.Value = 1;
            
            plotSettings = parsePlotSettings(app);
            
            [file,path] = uigetfile([app.outputDir.Value '\*.mat']);

            [~,~,~] = processWormdata(fullfile(path,file), plotSettings);

            app.validateRiseFall.Value = 1;


        end

        % Button pushed function: setTracePos
        function setTracePosButtonPushed(app, event)
            f = gcf; 
            app.TracesEditField.Value = num2str(f.Position);
        end

        % Button pushed function: setGraphPos
        function setGraphPosButtonPushed(app, event)
            f = gcf; 
            app.GraphsEditField.Value = num2str(f.Position);
        end

        % Button pushed function: PlotSingleTraceButton
        function PlotSingleTraceButtonPushed(app, event)
            prevDir = app.DataDir.Value;
            if isfolder(prevDir)
                searchStart = prevDir;
            else
                searchStart = app.defaultDataDir;
            end
 
            [file, path] = uigetfile([searchStart '\*.mat']);
            figure(app.IntestinalCalciumAppUIFigure)
            app.DataDir.Value = path;
            plotSettings = parsePlotSettings(app);
            plot_SingleTrace(fullfile(path,file), plotSettings)
        end

        % Menu selected function: saveSettings
        function saveSettingsMenuSelected(app, event)
            plotSettings = parsePlotSettings(app);
            uisave('plotSettings')
            figure(app.IntestinalCalciumAppUIFigure)
        end

        % Menu selected function: loadSettings
        function loadSettingsMenuSelected(app, event)
            uiload()
            figure(app.IntestinalCalciumAppUIFigure)
            applyPlotSettings(app,plotSettings);
        end

        % Menu selected function: savedefaultsettingsMenu
        function savedefaultsettingsMenuSelected(app, event)
            defaultPlotSettings = parsePlotSettings(app);
            [path, ~, ~] = fileparts(mfilename("fullpath"));
            save(fullfile(path, 'defaultPlotSettings.mat'),"defaultPlotSettings");

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create IntestinalCalciumAppUIFigure and hide until all components are created
            app.IntestinalCalciumAppUIFigure = uifigure('Visible', 'off');
            app.IntestinalCalciumAppUIFigure.Position = [100 100 473 475];
            app.IntestinalCalciumAppUIFigure.Name = 'Intestinal Calcium App';

            % Create FileMenu
            app.FileMenu = uimenu(app.IntestinalCalciumAppUIFigure);
            app.FileMenu.Text = 'File';

            % Create saveSettings
            app.saveSettings = uimenu(app.FileMenu);
            app.saveSettings.MenuSelectedFcn = createCallbackFcn(app, @saveSettingsMenuSelected, true);
            app.saveSettings.Text = 'save settings';

            % Create savedefaultsettingsMenu
            app.savedefaultsettingsMenu = uimenu(app.FileMenu);
            app.savedefaultsettingsMenu.MenuSelectedFcn = createCallbackFcn(app, @savedefaultsettingsMenuSelected, true);
            app.savedefaultsettingsMenu.Text = 'save default settings';

            % Create loadSettings
            app.loadSettings = uimenu(app.FileMenu);
            app.loadSettings.MenuSelectedFcn = createCallbackFcn(app, @loadSettingsMenuSelected, true);
            app.loadSettings.Separator = 'on';
            app.loadSettings.Text = 'load settings';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.IntestinalCalciumAppUIFigure);
            app.TabGroup.Position = [2 11 469 461];

            % Create TrackingTab
            app.TrackingTab = uitab(app.TabGroup);
            app.TrackingTab.Title = 'Tracking';

            % Create setTiffDirButton
            app.setTiffDirButton = uibutton(app.TrackingTab, 'push');
            app.setTiffDirButton.ButtonPushedFcn = createCallbackFcn(app, @setTiffDirButtonButtonPushed, true);
            app.setTiffDirButton.Position = [13 392 97 30];
            app.setTiffDirButton.Text = 'Set Tiff Dir';

            % Create tiffDir
            app.tiffDir = uieditfield(app.TrackingTab, 'text');
            app.tiffDir.HorizontalAlignment = 'right';
            app.tiffDir.FontSize = 10;
            app.tiffDir.Position = [119 394 252 28];
            app.tiffDir.Value = 'Select Tiff Directory';

            % Create IsRemoteCheckBox
            app.IsRemoteCheckBox = uicheckbox(app.TrackingTab);
            app.IsRemoteCheckBox.Tooltip = {'Check this box if data is on a server. This will temporarily copy Tif to local drive to analyze '};
            app.IsRemoteCheckBox.Text = 'Is Remote?';
            app.IsRemoteCheckBox.WordWrap = 'on';
            app.IsRemoteCheckBox.Position = [381 396 82 22];

            % Create LoadTiffCheckBox
            app.LoadTiffCheckBox = uicheckbox(app.TrackingTab);
            app.LoadTiffCheckBox.Tooltip = {'When checked, loads tiff into memory which is faster overall but requires sufficient RAM'};
            app.LoadTiffCheckBox.Text = 'Load Tiff';
            app.LoadTiffCheckBox.Position = [381 373 68 22];
            app.LoadTiffCheckBox.Value = true;

            % Create setRemoteDirButton
            app.setRemoteDirButton = uibutton(app.TrackingTab, 'push');
            app.setRemoteDirButton.ButtonPushedFcn = createCallbackFcn(app, @setRemoteDirButtonButtonPushed, true);
            app.setRemoteDirButton.Position = [13 338 97 30];
            app.setRemoteDirButton.Text = 'Set Remote Dir';

            % Create remoteDir
            app.remoteDir = uieditfield(app.TrackingTab, 'text');
            app.remoteDir.HorizontalAlignment = 'right';
            app.remoteDir.FontSize = 10;
            app.remoteDir.Position = [119 338 252 28];
            app.remoteDir.Value = 'Select Remote Directory';

            % Create UploadResultsCheckBox
            app.UploadResultsCheckBox = uicheckbox(app.TrackingTab);
            app.UploadResultsCheckBox.Text = 'Upload Results';
            app.UploadResultsCheckBox.WordWrap = 'on';
            app.UploadResultsCheckBox.Position = [380 335 87 30];
            app.UploadResultsCheckBox.Value = true;

            % Create DisplayOptionsPanel
            app.DisplayOptionsPanel = uipanel(app.TrackingTab);
            app.DisplayOptionsPanel.Title = 'Display Options';
            app.DisplayOptionsPanel.Position = [18 122 126 202];

            % Create PlotResultsCheckBox
            app.PlotResultsCheckBox = uicheckbox(app.DisplayOptionsPanel);
            app.PlotResultsCheckBox.Tooltip = {'Display the analysis progress in a figure. '};
            app.PlotResultsCheckBox.Text = 'Plot Results';
            app.PlotResultsCheckBox.Position = [3 158 86 22];
            app.PlotResultsCheckBox.Value = true;

            % Create ShowAxialSignalCheckBox
            app.ShowAxialSignalCheckBox = uicheckbox(app.DisplayOptionsPanel);
            app.ShowAxialSignalCheckBox.Tooltip = {'Display the analysis progress in a figure. '};
            app.ShowAxialSignalCheckBox.Text = 'Show Axial Signal';
            app.ShowAxialSignalCheckBox.Position = [3 137 118 22];
            app.ShowAxialSignalCheckBox.Value = true;

            % Create ShowNormalsCheckBox
            app.ShowNormalsCheckBox = uicheckbox(app.DisplayOptionsPanel);
            app.ShowNormalsCheckBox.Tooltip = {'Display the analysis progress in a figure. '};
            app.ShowNormalsCheckBox.Text = 'Show Normals';
            app.ShowNormalsCheckBox.Position = [3 113 100 22];
            app.ShowNormalsCheckBox.Value = true;

            % Create TroubleshootCheckBox
            app.TroubleshootCheckBox = uicheckbox(app.DisplayOptionsPanel);
            app.TroubleshootCheckBox.Tooltip = {'Show thresholding and results of morphological operations to troubleshoot tracking. '};
            app.TroubleshootCheckBox.Text = 'Troubleshoot';
            app.TroubleshootCheckBox.Position = [3 90 91 22];

            % Create RecordVideoCheckBox
            app.RecordVideoCheckBox = uicheckbox(app.DisplayOptionsPanel);
            app.RecordVideoCheckBox.Text = 'Record Video?';
            app.RecordVideoCheckBox.Position = [3 42 101 22];
            app.RecordVideoCheckBox.Value = true;

            % Create VideoFramerateEditField
            app.VideoFramerateEditField = uieditfield(app.DisplayOptionsPanel, 'numeric');
            app.VideoFramerateEditField.Tooltip = {'temporal downsample the analysis video.'; 'Record a frame every Nth frame of input tiff. '};
            app.VideoFramerateEditField.Position = [3 14 25 26];
            app.VideoFramerateEditField.Value = 15;

            % Create VideoFramerateEditFieldLabel
            app.VideoFramerateEditFieldLabel = uilabel(app.DisplayOptionsPanel);
            app.VideoFramerateEditFieldLabel.HorizontalAlignment = 'right';
            app.VideoFramerateEditFieldLabel.Position = [27 16 94 22];
            app.VideoFramerateEditFieldLabel.Text = 'Video Framerate';

            % Create ThresholdingButtonGroup
            app.ThresholdingButtonGroup = uibuttongroup(app.TrackingTab);
            app.ThresholdingButtonGroup.Title = 'Thresholding';
            app.ThresholdingButtonGroup.Position = [158 239 136 74];

            % Create AutoThresholdButton
            app.AutoThresholdButton = uiradiobutton(app.ThresholdingButtonGroup);
            app.AutoThresholdButton.Text = 'Auto Threshold';
            app.AutoThresholdButton.Position = [11 26 103 22];
            app.AutoThresholdButton.Value = true;

            % Create AdaptiveThresholdButton
            app.AdaptiveThresholdButton = uiradiobutton(app.ThresholdingButtonGroup);
            app.AdaptiveThresholdButton.Text = 'Adaptive Threshold';
            app.AdaptiveThresholdButton.Position = [11 4 125 22];

            % Create CropPixelsEditField
            app.CropPixelsEditField = uieditfield(app.TrackingTab, 'numeric');
            app.CropPixelsEditField.Tooltip = {'Number of edge pixels to crop. This helps remove black edges from split channel recordings. '};
            app.CropPixelsEditField.Position = [158 193 25 26];
            app.CropPixelsEditField.Value = 3;

            % Create CropPixelsEditFieldLabel
            app.CropPixelsEditFieldLabel = uilabel(app.TrackingTab);
            app.CropPixelsEditFieldLabel.HorizontalAlignment = 'right';
            app.CropPixelsEditFieldLabel.Position = [180 195 67 22];
            app.CropPixelsEditFieldLabel.Text = 'Crop Pixels';

            % Create InputFramerateEditField
            app.InputFramerateEditField = uieditfield(app.TrackingTab, 'numeric');
            app.InputFramerateEditField.Tooltip = {'Frame rate of the input tiff (frames per second)'};
            app.InputFramerateEditField.Position = [158 160 25 26];
            app.InputFramerateEditField.Value = 15;

            % Create InputFramerateEditFieldLabel
            app.InputFramerateEditFieldLabel = uilabel(app.TrackingTab);
            app.InputFramerateEditFieldLabel.HorizontalAlignment = 'right';
            app.InputFramerateEditFieldLabel.Position = [181 162 91 22];
            app.InputFramerateEditFieldLabel.Text = 'Input Framerate';

            % Create FlatfieldCorrectionEditField
            app.FlatfieldCorrectionEditField = uieditfield(app.TrackingTab, 'numeric');
            app.FlatfieldCorrectionEditField.Tooltip = {'Size of rolling ball used for flatfield correction (pixels)'};
            app.FlatfieldCorrectionEditField.Position = [158 124 25 26];
            app.FlatfieldCorrectionEditField.Value = 30;

            % Create FlatfieldCorrectionEditFieldLabel
            app.FlatfieldCorrectionEditFieldLabel = uilabel(app.TrackingTab);
            app.FlatfieldCorrectionEditFieldLabel.HorizontalAlignment = 'right';
            app.FlatfieldCorrectionEditFieldLabel.Position = [180 126 106 22];
            app.FlatfieldCorrectionEditFieldLabel.Text = 'Flatfield Correction';

            % Create StartFileSpinnerLabel
            app.StartFileSpinnerLabel = uilabel(app.TrackingTab);
            app.StartFileSpinnerLabel.HorizontalAlignment = 'right';
            app.StartFileSpinnerLabel.Position = [30 71 53 22];
            app.StartFileSpinnerLabel.Text = 'Start File';

            % Create StartFileSpinner
            app.StartFileSpinner = uispinner(app.TrackingTab);
            app.StartFileSpinner.Tooltip = {'Which tiff file in tiff Dir to begin analysis. '};
            app.StartFileSpinner.Position = [89 64 54 35];
            app.StartFileSpinner.Value = 1;

            % Create RunAnalysisButton
            app.RunAnalysisButton = uibutton(app.TrackingTab, 'push');
            app.RunAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @RunAnalysisButtonPushed, true);
            app.RunAnalysisButton.FontWeight = 'bold';
            app.RunAnalysisButton.Position = [156 17 271 82];
            app.RunAnalysisButton.Text = 'Run Analysis';

            % Create StartFrameSpinnerLabel
            app.StartFrameSpinnerLabel = uilabel(app.TrackingTab);
            app.StartFrameSpinnerLabel.HorizontalAlignment = 'right';
            app.StartFrameSpinnerLabel.Position = [17 25 68 22];
            app.StartFrameSpinnerLabel.Text = 'Start Frame';

            % Create StartFrameSpinner
            app.StartFrameSpinner = uispinner(app.TrackingTab);
            app.StartFrameSpinner.Tooltip = {'which frame to start the analysis from.'};
            app.StartFrameSpinner.Position = [89 18 54 35];
            app.StartFrameSpinner.Value = 1;

            % Create AnalyzeOASdataCheckBox
            app.AnalyzeOASdataCheckBox = uicheckbox(app.TrackingTab);
            app.AnalyzeOASdataCheckBox.ValueChangedFcn = createCallbackFcn(app, @AnalyzeOASdataCheckBox_2ValueChanged, true);
            app.AnalyzeOASdataCheckBox.Tooltip = {'Check this box if you are analyzing OAS data to run the proper tracking function'};
            app.AnalyzeOASdataCheckBox.Text = 'Analyze OAS data?';
            app.AnalyzeOASdataCheckBox.WordWrap = 'on';
            app.AnalyzeOASdataCheckBox.Position = [303 116 144 44];

            % Create PlottingTab
            app.PlottingTab = uitab(app.TabGroup);
            app.PlottingTab.Title = 'Plotting';

            % Create setDataDir
            app.setDataDir = uibutton(app.PlottingTab, 'push');
            app.setDataDir.ButtonPushedFcn = createCallbackFcn(app, @setDataDirButtonPushed, true);
            app.setDataDir.Position = [12 392 99 30];
            app.setDataDir.Text = 'Set Data Dir';

            % Create DataDir
            app.DataDir = uieditfield(app.PlottingTab, 'text');
            app.DataDir.HorizontalAlignment = 'right';
            app.DataDir.FontSize = 10;
            app.DataDir.Position = [119 394 252 28];
            app.DataDir.Value = 'Select directory containing wormdata.mat files';

            % Create outputDir
            app.outputDir = uieditfield(app.PlottingTab, 'text');
            app.outputDir.HorizontalAlignment = 'right';
            app.outputDir.FontSize = 10;
            app.outputDir.Position = [120 350 252 28];
            app.outputDir.Value = 'Select output directory if different from data dir';

            % Create setOutputDir
            app.setOutputDir = uibutton(app.PlottingTab, 'push');
            app.setOutputDir.ButtonPushedFcn = createCallbackFcn(app, @setOutputDirButtonPushed, true);
            app.setOutputDir.Position = [13 348 99 30];
            app.setOutputDir.Text = 'Set Output Dir';

            % Create controlnameEditFieldLabel
            app.controlnameEditFieldLabel = uilabel(app.PlottingTab);
            app.controlnameEditFieldLabel.HorizontalAlignment = 'right';
            app.controlnameEditFieldLabel.Position = [153 220 75 22];
            app.controlnameEditFieldLabel.Text = 'control name';

            % Create controlnameEditField
            app.controlnameEditField = uieditfield(app.PlottingTab, 'text');
            app.controlnameEditField.Tooltip = {'Specify what wormdata files to combine with mutant data as matched controls. this will create a mergedData.mat file containing mutant and control datasets for plotting.'};
            app.controlnameEditField.Position = [243 217 65 28];
            app.controlnameEditField.Value = 'wildtype';

            % Create PlotMatchedControlButton
            app.PlotMatchedControlButton = uibutton(app.PlottingTab, 'push');
            app.PlotMatchedControlButton.ButtonPushedFcn = createCallbackFcn(app, @PlotMatchedControlButtonPushed, true);
            app.PlotMatchedControlButton.FontWeight = 'bold';
            app.PlotMatchedControlButton.Position = [27 129 187 32];
            app.PlotMatchedControlButton.Text = 'Plot Matched Control';

            % Create MergeControlButton
            app.MergeControlButton = uibutton(app.PlottingTab, 'push');
            app.MergeControlButton.ButtonPushedFcn = createCallbackFcn(app, @MergeControlButtonPushed, true);
            app.MergeControlButton.FontWeight = 'bold';
            app.MergeControlButton.Position = [135 252 187 32];
            app.MergeControlButton.Text = 'Merge Control';

            % Create CombinewormDataButton
            app.CombinewormDataButton = uibutton(app.PlottingTab, 'push');
            app.CombinewormDataButton.ButtonPushedFcn = createCallbackFcn(app, @CombinewormDataButtonPushed, true);
            app.CombinewormDataButton.FontWeight = 'bold';
            app.CombinewormDataButton.Tooltip = {'This will seach through folders '};
            app.CombinewormDataButton.Position = [135 298 187 32];
            app.CombinewormDataButton.Text = 'Combine wormData';

            % Create Panel
            app.Panel = uipanel(app.PlottingTab);
            app.Panel.Position = [16 11 212 102];

            % Create PlotSingleTraceButton
            app.PlotSingleTraceButton = uibutton(app.Panel, 'push');
            app.PlotSingleTraceButton.ButtonPushedFcn = createCallbackFcn(app, @PlotSingleTraceButtonPushed, true);
            app.PlotSingleTraceButton.FontWeight = 'bold';
            app.PlotSingleTraceButton.Position = [10 50 187 32];
            app.PlotSingleTraceButton.Text = 'Plot Single Trace';

            % Create spikeEditFieldLabel
            app.spikeEditFieldLabel = uilabel(app.Panel);
            app.spikeEditFieldLabel.HorizontalAlignment = 'right';
            app.spikeEditFieldLabel.Position = [12 14 43 22];
            app.spikeEditFieldLabel.Text = 'spike #';

            % Create spikeEditField
            app.spikeEditField = uieditfield(app.Panel, 'numeric');
            app.spikeEditField.Position = [66 14 27 22];
            app.spikeEditField.Value = 1;

            % Create windowsizeLabel
            app.windowsizeLabel = uilabel(app.Panel);
            app.windowsizeLabel.HorizontalAlignment = 'center';
            app.windowsizeLabel.Position = [110 10 48 30];
            app.windowsizeLabel.Text = {'window '; 'size'};

            % Create windowsizeEditField
            app.windowsizeEditField = uieditfield(app.Panel, 'numeric');
            app.windowsizeEditField.Position = [166 14 27 22];
            app.windowsizeEditField.Value = 10;

            % Create Panel_2
            app.Panel_2 = uipanel(app.PlottingTab);
            app.Panel_2.Position = [238 11 222 156];

            % Create PlotMultipleGenotypesButton
            app.PlotMultipleGenotypesButton = uibutton(app.Panel_2, 'push');
            app.PlotMultipleGenotypesButton.ButtonPushedFcn = createCallbackFcn(app, @PlotMultipleGenotypesButtonPushed, true);
            app.PlotMultipleGenotypesButton.FontWeight = 'bold';
            app.PlotMultipleGenotypesButton.Position = [21 114 187 32];
            app.PlotMultipleGenotypesButton.Text = 'Plot Multiple Genotypes';

            % Create mutantgenotypesEditFieldLabel
            app.mutantgenotypesEditFieldLabel = uilabel(app.Panel_2);
            app.mutantgenotypesEditFieldLabel.HorizontalAlignment = 'right';
            app.mutantgenotypesEditFieldLabel.Position = [65 90 100 22];
            app.mutantgenotypesEditFieldLabel.Text = 'mutant genotypes';

            % Create mutantgenotypesEditField
            app.mutantgenotypesEditField = uieditfield(app.Panel_2, 'text');
            app.mutantgenotypesEditField.Tooltip = {'list genotypes to plot together. Must match the filenames of merged data files (eg wildtype would correspond to the file "wildtype_mergedData.mat")'};
            app.mutantgenotypesEditField.Position = [23 65 183 28];

            % Create ColumnsEditFieldLabel
            app.ColumnsEditFieldLabel = uilabel(app.Panel_2);
            app.ColumnsEditFieldLabel.HorizontalAlignment = 'right';
            app.ColumnsEditFieldLabel.Position = [62 38 63 22];
            app.ColumnsEditFieldLabel.Text = '# Columns';

            % Create numColumns
            app.numColumns = uieditfield(app.Panel_2, 'numeric');
            app.numColumns.Tooltip = {'Hint - set to zero to have the same number of columns as genotypes. '};
            app.numColumns.Position = [137 38 31 22];

            % Create PrefixEditFieldLabel
            app.PrefixEditFieldLabel = uilabel(app.Panel_2);
            app.PrefixEditFieldLabel.HorizontalAlignment = 'right';
            app.PrefixEditFieldLabel.Position = [49 9 39 22];
            app.PrefixEditFieldLabel.Text = 'Prefix ';

            % Create PrefixEditField
            app.PrefixEditField = uieditfield(app.Panel_2, 'text');
            app.PrefixEditField.Position = [103 6 83 28];
            app.PrefixEditField.Value = 'Mutants';

            % Create PlotSettingsTab
            app.PlotSettingsTab = uitab(app.TabGroup);
            app.PlotSettingsTab.Title = 'Plot Settings';

            % Create MultigenotypefigurepositionPanel
            app.MultigenotypefigurepositionPanel = uipanel(app.PlotSettingsTab);
            app.MultigenotypefigurepositionPanel.TitlePosition = 'centertop';
            app.MultigenotypefigurepositionPanel.Title = 'Multi-genotype figure position';
            app.MultigenotypefigurepositionPanel.FontWeight = 'bold';
            app.MultigenotypefigurepositionPanel.Position = [203 9 251 79];

            % Create GraphsEditFieldLabel
            app.GraphsEditFieldLabel = uilabel(app.MultigenotypefigurepositionPanel);
            app.GraphsEditFieldLabel.HorizontalAlignment = 'right';
            app.GraphsEditFieldLabel.Position = [5 32 44 22];
            app.GraphsEditFieldLabel.Text = 'Graphs';

            % Create GraphsEditField
            app.GraphsEditField = uieditfield(app.MultigenotypefigurepositionPanel, 'text');
            app.GraphsEditField.Position = [64 34 106 18];
            app.GraphsEditField.Value = '253 581 600 156';

            % Create TracesEditFieldLabel
            app.TracesEditFieldLabel = uilabel(app.MultigenotypefigurepositionPanel);
            app.TracesEditFieldLabel.HorizontalAlignment = 'right';
            app.TracesEditFieldLabel.Position = [8 6 41 22];
            app.TracesEditFieldLabel.Text = 'Traces';

            % Create TracesEditField
            app.TracesEditField = uieditfield(app.MultigenotypefigurepositionPanel, 'text');
            app.TracesEditField.Position = [64 8 106 18];
            app.TracesEditField.Value = '147 281 1205 396';

            % Create setGraphPos
            app.setGraphPos = uibutton(app.MultigenotypefigurepositionPanel, 'push');
            app.setGraphPos.ButtonPushedFcn = createCallbackFcn(app, @setGraphPosButtonPushed, true);
            app.setGraphPos.Position = [190 32 46 22];
            app.setGraphPos.Text = 'Set';

            % Create setTracePos
            app.setTracePos = uibutton(app.MultigenotypefigurepositionPanel, 'push');
            app.setTracePos.ButtonPushedFcn = createCallbackFcn(app, @setTracePosButtonPushed, true);
            app.setTracePos.Position = [190 6 46 22];
            app.setTracePos.Text = 'Set';

            % Create YAxesPanel
            app.YAxesPanel = uipanel(app.PlotSettingsTab);
            app.YAxesPanel.Title = 'Y Axes';
            app.YAxesPanel.FontWeight = 'bold';
            app.YAxesPanel.Position = [24 335 174 87];

            % Create AxialSignallimLabel
            app.AxialSignallimLabel = uilabel(app.YAxesPanel);
            app.AxialSignallimLabel.HorizontalAlignment = 'right';
            app.AxialSignallimLabel.Position = [-2 41 96 22];
            app.AxialSignallimLabel.Text = 'Axial Signal Lims';

            % Create axialYLim
            app.axialYLim = uieditfield(app.YAxesPanel, 'text');
            app.axialYLim.HorizontalAlignment = 'right';
            app.axialYLim.Tooltip = {'Y Axis Limits [Min Max] (arbitrary units) for Axial signal plots.'};
            app.axialYLim.Position = [97 41 74 22];
            app.axialYLim.Value = '-500 12000';

            % Create BulkSignallimsLabel
            app.BulkSignallimsLabel = uilabel(app.YAxesPanel);
            app.BulkSignallimsLabel.HorizontalAlignment = 'right';
            app.BulkSignallimsLabel.Position = [-2 8 94 22];
            app.BulkSignallimsLabel.Text = 'Bulk Signal Lims';

            % Create bulkYLim
            app.bulkYLim = uieditfield(app.YAxesPanel, 'text');
            app.bulkYLim.HorizontalAlignment = 'right';
            app.bulkYLim.Tooltip = {'Y Axis Limits [Min Max] (arbitrary units) for Bulk Signal plots.'};
            app.bulkYLim.Position = [98 8 74 22];
            app.bulkYLim.Value = '-500 5000';

            % Create XAxesPanel
            app.XAxesPanel = uipanel(app.PlotSettingsTab);
            app.XAxesPanel.Title = 'X Axes';
            app.XAxesPanel.FontWeight = 'bold';
            app.XAxesPanel.Position = [24 239 174 84];

            % Create BulkSignallimitsLabel
            app.BulkSignallimitsLabel = uilabel(app.XAxesPanel);
            app.BulkSignallimitsLabel.HorizontalAlignment = 'right';
            app.BulkSignallimitsLabel.Position = [-2 36 94 22];
            app.BulkSignallimitsLabel.Text = 'Bulk Signal Lims';

            % Create bulkXLim
            app.bulkXLim = uieditfield(app.XAxesPanel, 'text');
            app.bulkXLim.HorizontalAlignment = 'right';
            app.bulkXLim.Tooltip = {'X axis limits (minutes) for bulk signal traces.'};
            app.bulkXLim.Position = [130 36 39 22];
            app.bulkXLim.Value = '0 10';

            % Create BulkSignallimitsLabel_2
            app.BulkSignallimitsLabel_2 = uilabel(app.XAxesPanel);
            app.BulkSignallimitsLabel_2.HorizontalAlignment = 'right';
            app.BulkSignallimitsLabel_2.Position = [-1 10 99 22];
            app.BulkSignallimitsLabel_2.Text = 'Axial Signal Ticks';

            % Create axialXTickInt
            app.axialXTickInt = uieditfield(app.XAxesPanel, 'text');
            app.axialXTickInt.HorizontalAlignment = 'right';
            app.axialXTickInt.Tooltip = {'Increments for X Axis Ticks (minutes) for Axial Signal plots'};
            app.axialXTickInt.Position = [130 10 39 22];
            app.axialXTickInt.Value = '2';

            % Create SortTracesByButtonGroup
            app.SortTracesByButtonGroup = uibuttongroup(app.PlotSettingsTab);
            app.SortTracesByButtonGroup.Title = 'Sort Traces By:';
            app.SortTracesByButtonGroup.FontWeight = 'bold';
            app.SortTracesByButtonGroup.Position = [207 280 120 141];

            % Create dontSort
            app.dontSort = uiradiobutton(app.SortTracesByButtonGroup);
            app.dontSort.Text = 'Dont sort';
            app.dontSort.Position = [13 101 71 22];
            app.dontSort.Value = true;

            % Create sortFreq
            app.sortFreq = uiradiobutton(app.SortTracesByButtonGroup);
            app.sortFreq.Text = 'Spike frequency';
            app.sortFreq.Position = [13 79 108 22];

            % Create sortAmp
            app.sortAmp = uiradiobutton(app.SortTracesByButtonGroup);
            app.sortAmp.Text = 'Spike amplitude';
            app.sortAmp.Position = [12 57 107 22];

            % Create HistogramBinsPanel
            app.HistogramBinsPanel = uipanel(app.PlotSettingsTab);
            app.HistogramBinsPanel.Title = 'Histogram Bins';
            app.HistogramBinsPanel.FontWeight = 'bold';
            app.HistogramBinsPanel.Position = [360 309 100 112];

            % Create MinEditFieldLabel
            app.MinEditFieldLabel = uilabel(app.HistogramBinsPanel);
            app.MinEditFieldLabel.HorizontalAlignment = 'right';
            app.MinEditFieldLabel.Position = [36 57 25 22];
            app.MinEditFieldLabel.Text = 'Min';

            % Create MinEditField
            app.MinEditField = uieditfield(app.HistogramBinsPanel, 'numeric');
            app.MinEditField.Position = [65 57 28 22];

            % Create MaxEditFieldLabel
            app.MaxEditFieldLabel = uilabel(app.HistogramBinsPanel);
            app.MaxEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxEditFieldLabel.Position = [35 31 27 22];
            app.MaxEditFieldLabel.Text = 'Max';

            % Create MaxEditField
            app.MaxEditField = uieditfield(app.HistogramBinsPanel, 'numeric');
            app.MaxEditField.Position = [66 31 28 22];
            app.MaxEditField.Value = 90;

            % Create IncrementEditFieldLabel
            app.IncrementEditFieldLabel = uilabel(app.HistogramBinsPanel);
            app.IncrementEditFieldLabel.HorizontalAlignment = 'right';
            app.IncrementEditFieldLabel.Position = [3 5 59 22];
            app.IncrementEditFieldLabel.Text = 'Increment';

            % Create IncrementEditField
            app.IncrementEditField = uieditfield(app.HistogramBinsPanel, 'numeric');
            app.IncrementEditField.Position = [66 5 28 22];
            app.IncrementEditField.Value = 2;

            % Create PeakDetectionPanel
            app.PeakDetectionPanel = uipanel(app.PlotSettingsTab);
            app.PeakDetectionPanel.Title = 'Peak Detection';
            app.PeakDetectionPanel.FontWeight = 'bold';
            app.PeakDetectionPanel.Position = [24 90 174 143];

            % Create PeakProminenceAULabel
            app.PeakProminenceAULabel = uilabel(app.PeakDetectionPanel);
            app.PeakProminenceAULabel.HorizontalAlignment = 'right';
            app.PeakProminenceAULabel.Position = [2 95 127 22];
            app.PeakProminenceAULabel.Text = 'Peak Prominence (AU)';

            % Create peakthreshold
            app.peakthreshold = uieditfield(app.PeakDetectionPanel, 'numeric');
            app.peakthreshold.Tooltip = {'Min Peak Prominence. See "findpeaks()".'};
            app.peakthreshold.Position = [132 95 37 22];
            app.peakthreshold.Value = 750;

            % Create MinDistanceframesEditFieldLabel
            app.MinDistanceframesEditFieldLabel = uilabel(app.PeakDetectionPanel);
            app.MinDistanceframesEditFieldLabel.HorizontalAlignment = 'right';
            app.MinDistanceframesEditFieldLabel.Position = [2 68 123 22];
            app.MinDistanceframesEditFieldLabel.Text = 'Min Distance (frames)';

            % Create peakdistance
            app.peakdistance = uieditfield(app.PeakDetectionPanel, 'numeric');
            app.peakdistance.Tooltip = {'Min Peak Distance. See "findpeaks()".'};
            app.peakdistance.Position = [132 68 37 22];
            app.peakdistance.Value = 150;

            % Create MinWidthframesEditFieldLabel
            app.MinWidthframesEditFieldLabel = uilabel(app.PeakDetectionPanel);
            app.MinWidthframesEditFieldLabel.HorizontalAlignment = 'right';
            app.MinWidthframesEditFieldLabel.Position = [3 41 107 22];
            app.MinWidthframesEditFieldLabel.Text = 'Min Width (frames)';

            % Create peakwidth
            app.peakwidth = uieditfield(app.PeakDetectionPanel, 'numeric');
            app.peakwidth.Tooltip = {'Min Peak Width. See "findpeaks()".'};
            app.peakwidth.Position = [132 41 37 22];
            app.peakwidth.Value = 15;

            % Create PlotwindowsecLabel
            app.PlotwindowsecLabel = uilabel(app.PeakDetectionPanel);
            app.PlotwindowsecLabel.HorizontalAlignment = 'right';
            app.PlotwindowsecLabel.Position = [1 12 102 22];
            app.PlotwindowsecLabel.Text = 'Plot Window (sec)';

            % Create spikeProfileWindow
            app.spikeProfileWindow = uieditfield(app.PeakDetectionPanel, 'numeric');
            app.spikeProfileWindow.Position = [130 12 37 22];
            app.spikeProfileWindow.Value = 30;

            % Create FrameRateEditFieldLabel
            app.FrameRateEditFieldLabel = uilabel(app.PlotSettingsTab);
            app.FrameRateEditFieldLabel.HorizontalAlignment = 'right';
            app.FrameRateEditFieldLabel.Position = [360 199 68 22];
            app.FrameRateEditFieldLabel.Text = 'Frame Rate';

            % Create FrameRateEditField
            app.FrameRateEditField = uieditfield(app.PlotSettingsTab, 'numeric');
            app.FrameRateEditField.Tooltip = {'Frames per second of recorded data'};
            app.FrameRateEditField.Position = [431 201 29 19];
            app.FrameRateEditField.Value = 15;

            % Create NormalizeCheckBox
            app.NormalizeCheckBox = uicheckbox(app.PlotSettingsTab);
            app.NormalizeCheckBox.ValueChangedFcn = createCallbackFcn(app, @NormalizeCheckBoxValueChanged, true);
            app.NormalizeCheckBox.Text = 'Normalize Bulk Signal';
            app.NormalizeCheckBox.Position = [207 235 154 22];

            % Create AutoFixAxialSignalCheckBox
            app.AutoFixAxialSignalCheckBox = uicheckbox(app.PlotSettingsTab);
            app.AutoFixAxialSignalCheckBox.Tooltip = {'Re-analyze axial signal to correct for head tail flips.'};
            app.AutoFixAxialSignalCheckBox.Text = 'Auto-Fix Axial Signal';
            app.AutoFixAxialSignalCheckBox.Position = [207 139 132 24];
            app.AutoFixAxialSignalCheckBox.Value = true;

            % Create toQuerry
            app.toQuerry = uieditfield(app.PlotSettingsTab, 'numeric');
            app.toQuerry.Tooltip = {'fraction of axial signal to compare for head/tail reassignment'};
            app.toQuerry.Position = [207 118 38 19];
            app.toQuerry.Value = 0.15;

            % Create FractiontoQuerrylengthLabel
            app.FractiontoQuerrylengthLabel = uilabel(app.PlotSettingsTab);
            app.FractiontoQuerrylengthLabel.HorizontalAlignment = 'center';
            app.FractiontoQuerrylengthLabel.WordWrap = 'on';
            app.FractiontoQuerrylengthLabel.Position = [242 118 118 19];
            app.FractiontoQuerrylengthLabel.Text = '% of length to querry';

            % Create SortDirectionButtonGroup
            app.SortDirectionButtonGroup = uibuttongroup(app.PlotSettingsTab);
            app.SortDirectionButtonGroup.Title = 'Sort Direction';
            app.SortDirectionButtonGroup.FontWeight = 'bold';
            app.SortDirectionButtonGroup.Position = [207 259 120 79];

            % Create descendButton
            app.descendButton = uiradiobutton(app.SortDirectionButtonGroup);
            app.descendButton.Text = 'descend';
            app.descendButton.Position = [10 38 67 22];
            app.descendButton.Value = true;

            % Create ascendButton
            app.ascendButton = uiradiobutton(app.SortDirectionButtonGroup);
            app.ascendButton.Text = 'ascend';
            app.ascendButton.Position = [10 20 65 22];

            % Create shuffleButton
            app.shuffleButton = uiradiobutton(app.SortDirectionButtonGroup);
            app.shuffleButton.Tooltip = {'Only works with "dont sort"'};
            app.shuffleButton.Text = 'shuffle';
            app.shuffleButton.Position = [10 1 65 22];

            % Create BulkSignalPlotsButtonGroup
            app.BulkSignalPlotsButtonGroup = uibuttongroup(app.PlotSettingsTab);
            app.BulkSignalPlotsButtonGroup.Tooltip = {'should bulk signal be plotted in its own axis or overlaid on the axialSignal'};
            app.BulkSignalPlotsButtonGroup.Title = 'Bulk Signal Plots';
            app.BulkSignalPlotsButtonGroup.FontWeight = 'bold';
            app.BulkSignalPlotsButtonGroup.FontSize = 11;
            app.BulkSignalPlotsButtonGroup.Position = [361 231 100 70];

            % Create SeparateButton
            app.SeparateButton = uiradiobutton(app.BulkSignalPlotsButtonGroup);
            app.SeparateButton.Text = 'Separate';
            app.SeparateButton.Position = [11 25 71 22];
            app.SeparateButton.Value = true;

            % Create OverlayButton
            app.OverlayButton = uiradiobutton(app.BulkSignalPlotsButtonGroup);
            app.OverlayButton.Text = 'Overlay';
            app.OverlayButton.Position = [11 3 65 22];

            % Create PlotsEditFieldLabel
            app.PlotsEditFieldLabel = uilabel(app.PlotSettingsTab);
            app.PlotsEditFieldLabel.HorizontalAlignment = 'right';
            app.PlotsEditFieldLabel.Position = [383 166 42 22];
            app.PlotsEditFieldLabel.Text = '# Plots';

            % Create numPlotsEditField
            app.numPlotsEditField = uieditfield(app.PlotSettingsTab, 'numeric');
            app.numPlotsEditField.Tooltip = {'number of bulk/axial signal plots to display. does not affect other graphs. '};
            app.numPlotsEditField.Position = [431 166 29 22];
            app.numPlotsEditField.Value = 7;

            % Create PlotMatchedControlButton_2
            app.PlotMatchedControlButton_2 = uibutton(app.PlotSettingsTab, 'push');
            app.PlotMatchedControlButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotMatchedControlButtonPushed, true);
            app.PlotMatchedControlButton_2.FontWeight = 'bold';
            app.PlotMatchedControlButton_2.Position = [30 46 162 32];
            app.PlotMatchedControlButton_2.Text = 'Plot Matched Control';

            % Create PlotMultipleGenotypesButton_2
            app.PlotMultipleGenotypesButton_2 = uibutton(app.PlotSettingsTab, 'push');
            app.PlotMultipleGenotypesButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotMultipleGenotypesButtonPushed, true);
            app.PlotMultipleGenotypesButton_2.FontWeight = 'bold';
            app.PlotMultipleGenotypesButton_2.Position = [30 8 162 32];
            app.PlotMultipleGenotypesButton_2.Text = 'Plot Multiple Genotypes';

            % Create AnalyzeOASdataCheckBox_2
            app.AnalyzeOASdataCheckBox_2 = uicheckbox(app.PlotSettingsTab);
            app.AnalyzeOASdataCheckBox_2.ValueChangedFcn = createCallbackFcn(app, @AnalyzeOASdataCheckBox_2ValueChanged, true);
            app.AnalyzeOASdataCheckBox_2.Tooltip = {'Check if working with OpenAutoScope data. This will toggle several required settings'};
            app.AnalyzeOASdataCheckBox_2.Text = 'Analyze OAS data?';
            app.AnalyzeOASdataCheckBox_2.WordWrap = 'on';
            app.AnalyzeOASdataCheckBox_2.Position = [207 207 144 18];

            % Create EqualizeExpDurationCheckBox
            app.EqualizeExpDurationCheckBox = uicheckbox(app.PlotSettingsTab);
            app.EqualizeExpDurationCheckBox.Tooltip = {'check to trim all experiments to the length of the shortest recording. Required for OAS.'};
            app.EqualizeExpDurationCheckBox.Text = 'Equalize Exp Duration?';
            app.EqualizeExpDurationCheckBox.WordWrap = 'on';
            app.EqualizeExpDurationCheckBox.Position = [207 188 144 14];

            % Create AxialSignalColorLabel
            app.AxialSignalColorLabel = uilabel(app.PlotSettingsTab);
            app.AxialSignalColorLabel.HorizontalAlignment = 'center';
            app.AxialSignalColorLabel.Position = [381 129 71 30];
            app.AxialSignalColorLabel.Text = {'Axial Signal '; 'Color'};

            % Create AxialSignalColorDropDown
            app.AxialSignalColorDropDown = uidropdown(app.PlotSettingsTab);
            app.AxialSignalColorDropDown.Items = {'viridis', 'turbo', 'bone', 'jet', 'inferno', 'magma', 'plasma', 'twilight', 'cividis'};
            app.AxialSignalColorDropDown.Position = [376 99 81 28];
            app.AxialSignalColorDropDown.Value = 'viridis';

            % Create miscsettingsTab
            app.miscsettingsTab = uitab(app.TabGroup);
            app.miscsettingsTab.Title = 'misc settings';

            % Create ColorsPanel
            app.ColorsPanel = uipanel(app.miscsettingsTab);
            app.ColorsPanel.Title = 'Colors';
            app.ColorsPanel.FontWeight = 'bold';
            app.ColorsPanel.Position = [216 306 233 115];

            % Create mutantcolorEditFieldLabel
            app.mutantcolorEditFieldLabel = uilabel(app.ColorsPanel);
            app.mutantcolorEditFieldLabel.HorizontalAlignment = 'right';
            app.mutantcolorEditFieldLabel.Position = [43 67 71 22];
            app.mutantcolorEditFieldLabel.Text = 'mutant color';

            % Create mtColor
            app.mtColor = uieditfield(app.ColorsPanel, 'text');
            app.mtColor.Position = [125 67 100 22];
            app.mtColor.Value = '0.66 0.74 0.91';

            % Create wildtypecolorEditFieldLabel
            app.wildtypecolorEditFieldLabel = uilabel(app.ColorsPanel);
            app.wildtypecolorEditFieldLabel.HorizontalAlignment = 'right';
            app.wildtypecolorEditFieldLabel.Position = [34 36 81 22];
            app.wildtypecolorEditFieldLabel.Text = 'wild type color';

            % Create wtEdgeColor
            app.wtEdgeColor = uieditfield(app.ColorsPanel, 'text');
            app.wtEdgeColor.Position = [125 36 100 22];
            app.wtEdgeColor.Value = '0.4 0.4 0.4';

            % Create wildtypeedgecolorEditFieldLabel
            app.wildtypeedgecolorEditFieldLabel = uilabel(app.ColorsPanel);
            app.wildtypeedgecolorEditFieldLabel.HorizontalAlignment = 'center';
            app.wildtypeedgecolorEditFieldLabel.WordWrap = 'on';
            app.wildtypeedgecolorEditFieldLabel.Position = [6 3 114 30];
            app.wildtypeedgecolorEditFieldLabel.Text = 'wild type edge color';

            % Create wtColor
            app.wtColor = uieditfield(app.ColorsPanel, 'text');
            app.wtColor.Position = [126 7 100 22];
            app.wtColor.Value = '0.26 0.34 0.51';

            % Create validateRiseFallButton
            app.validateRiseFallButton = uibutton(app.miscsettingsTab, 'push');
            app.validateRiseFallButton.ButtonPushedFcn = createCallbackFcn(app, @validateRiseFallButtonPushed, true);
            app.validateRiseFallButton.Position = [67 291 104 36];
            app.validateRiseFallButton.Text = 'Validate Rise/Fall';

            % Create SpikeKineticsPanel
            app.SpikeKineticsPanel = uipanel(app.miscsettingsTab);
            app.SpikeKineticsPanel.Title = 'Spike Kinetics ';
            app.SpikeKineticsPanel.FontWeight = 'bold';
            app.SpikeKineticsPanel.Position = [50 334 143 87];

            % Create showFitParams
            app.showFitParams = uicheckbox(app.SpikeKineticsPanel);
            app.showFitParams.Text = 'Show rise/fall time fit?';
            app.showFitParams.Position = [3 44 139 22];

            % Create validatePropagationRate
            app.validatePropagationRate = uicheckbox(app.SpikeKineticsPanel);
            app.validatePropagationRate.Text = 'Validate propagation?';
            app.validatePropagationRate.Position = [1 2 138 24];

            % Create validateRiseFall
            app.validateRiseFall = uicheckbox(app.SpikeKineticsPanel);
            app.validateRiseFall.Text = 'Validate rise/fall?';
            app.validateRiseFall.Position = [2 25 138 24];

            % Show the figure after all components are created
            app.IntestinalCalciumAppUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Intestinal_Calcium_app_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.IntestinalCalciumAppUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.IntestinalCalciumAppUIFigure)
        end
    end
end