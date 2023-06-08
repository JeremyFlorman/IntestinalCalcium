classdef Intestinal_Calcium_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        IntestinalCalciumAppUIFigure   matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        TrackingTab                    matlab.ui.container.Tab
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
        numColumns                     matlab.ui.control.NumericEditField
        ColumnsEditFieldLabel          matlab.ui.control.Label
        PlotMultipleGenotypesButton    matlab.ui.control.Button
        CombinewormDataButton          matlab.ui.control.Button
        MergeControlButton             matlab.ui.control.Button
        PlotMatchedControlButton       matlab.ui.control.Button
        controlnameEditField           matlab.ui.control.EditField
        controlnameEditFieldLabel      matlab.ui.control.Label
        mutantgenotypesEditField       matlab.ui.control.EditField
        mutantgenotypesEditFieldLabel  matlab.ui.control.Label
        setOutputDir                   matlab.ui.control.Button
        outputDir                      matlab.ui.control.EditField
        DataDir                        matlab.ui.control.EditField
        setDataDir                     matlab.ui.control.Button
        PlotSettingsTab                matlab.ui.container.Tab
        numPlotsEditField              matlab.ui.control.NumericEditField
        PlotsEditFieldLabel            matlab.ui.control.Label
        BulkSignalPlotsButtonGroup     matlab.ui.container.ButtonGroup
        OverlayButton                  matlab.ui.control.RadioButton
        SeparateButton                 matlab.ui.control.RadioButton
        SortDirectionButtonGroup       matlab.ui.container.ButtonGroup
        ascendButton                   matlab.ui.control.RadioButton
        descendButton                  matlab.ui.control.RadioButton
        FractiontoQuerrylengthLabel    matlab.ui.control.Label
        toQuerry                       matlab.ui.control.NumericEditField
        test                           matlab.ui.control.Button
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
        ColorsPanel                    matlab.ui.container.Panel
        wtColor                        matlab.ui.control.EditField
        wildtypeedgecolorEditFieldLabel  matlab.ui.control.Label
        wtEdgeColor                    matlab.ui.control.EditField
        wildtypecolorEditFieldLabel    matlab.ui.control.Label
        mtColor                        matlab.ui.control.EditField
        mutantcolorEditFieldLabel      matlab.ui.control.Label
        YAxesPanel                     matlab.ui.container.Panel
        bulkYLim                       matlab.ui.control.EditField
        BulkSignallimsLabel            matlab.ui.control.Label
        axialYLim                      matlab.ui.control.EditField
        AxialSignallimLabel            matlab.ui.control.Label
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


            plotSettings.wtcolor = str2num(app.wtColor.Value);
            plotSettings.mtcolor = str2num(app.mtColor.Value);
            plotSettings.mtedgecolor = str2num(app.wtEdgeColor.Value);



            plotSettings.binedges = app.MinEditField.Value:app.IncrementEditField.Value:app.MaxEditField.Value;


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
            end

            if app.SeparateButton.Value
                plotSettings.overlayplots = 0;
            elseif app.OverlayButton.Value
                plotSettings.overlayplots = 1;
            end

            plotSettings.tolimit = app.numPlotsEditField.Value;
            plotSettings.controlname = app.controlnameEditField.Value;
            
            
            gtyps = split(app.mutantgenotypesEditField.Value, {' ' , ',', ';'});
            gtyps = gtyps(~cellfun(@isempty,gtyps));

            if length(gtyps) == 1
                gtyps = gtyps{:};
            end

            plotSettings.genotypes = gtyps;


            if app.numColumns.Value == 0 
                numcolumns = length(gtyps);
            else
                numcolumns = app.numColumns.Value;
            end

            
            plotSettings.numColumns = numcolumns;

        end




    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: setTiffDirButton
        function setTiffDirButtonButtonPushed(app, event)
            prevDir = app.tiffDir.Value;
            if isfolder(prevDir)
                searchStart = prevDir;
            else
                searchStart = 'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants';
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
                searchStart = 'Y:\';
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
            freelyMovingAnalysis_Func(parsedInputs)

        end

        % Button pushed function: setDataDir
        function setDataDirButtonPushed(app, event)
            prevDir = app.DataDir.Value;
            if isfolder(prevDir)
                searchStart = prevDir;
            else
                searchStart = 'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants';
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
                searchStart = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\';
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
                datadir = uigetdir('Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants',...
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
                    outputdir = uigetdir('Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants',...
                        'Oops! looks like you forgot to select your data directory!');
                    app.DataDir.Value = outputdir;
                    figure(app.IntestinalCalciumAppUIFigure)
                end
            end
            mergeControl(outputdir, app.controlnameEditField.Value)
        end

        % Button pushed function: test
        function testButtonPushed(app, event)
            settings = parsePlotSettings(app);
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

        % Button pushed function: PlotMatchedControlButton
        function PlotMatchedControlButtonPushed(app, event)
            outputdir = app.outputDir.Value;
            if ~isfolder(outputdir)
                if isfolder(app.DataDir.Value)
                    outputdir = app.DataDir.Value;
                else
                    outputdir = uigetdir('Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants',...
                        'Oops! looks like you forgot to select your data directory!');
                    app.DataDir.Value = outputdir;
                    figure(app.IntestinalCalciumAppUIFigure)
                end
            end
            settings = parsePlotSettings(app);
            plot_MatchedControl(app.outputDir.Value,settings);

        end

        % Button pushed function: PlotMultipleGenotypesButton
        function PlotMultipleGenotypesButtonPushed(app, event)
            
            outputdir = app.outputDir.Value;
            if ~isfolder(outputdir)
                if isfolder(app.DataDir.Value)
                    outputdir = app.DataDir.Value;
                else
                    outputdir = uigetdir('Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants',...
                        'Oops! looks like you forgot to select your data directory!');
                    app.DataDir.Value = outputdir;
                    figure(app.IntestinalCalciumAppUIFigure)
                end
            end
            settings = parsePlotSettings(app);
            plot_MultiGenotype(outputdir, settings)


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

            % Create mutantgenotypesEditFieldLabel
            app.mutantgenotypesEditFieldLabel = uilabel(app.PlottingTab);
            app.mutantgenotypesEditFieldLabel.HorizontalAlignment = 'right';
            app.mutantgenotypesEditFieldLabel.Position = [177 75 100 22];
            app.mutantgenotypesEditFieldLabel.Text = 'mutant genotypes';

            % Create mutantgenotypesEditField
            app.mutantgenotypesEditField = uieditfield(app.PlottingTab, 'text');
            app.mutantgenotypesEditField.Tooltip = {'list genotypes to plot together. Must match the filenames of merged data files (eg wildtype would correspond to the file "wildtype_mergedData.mat")'};
            app.mutantgenotypesEditField.Position = [135 50 183 28];

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
            app.PlotMatchedControlButton.Position = [133 158 187 32];
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

            % Create PlotMultipleGenotypesButton
            app.PlotMultipleGenotypesButton = uibutton(app.PlottingTab, 'push');
            app.PlotMultipleGenotypesButton.ButtonPushedFcn = createCallbackFcn(app, @PlotMultipleGenotypesButtonPushed, true);
            app.PlotMultipleGenotypesButton.FontWeight = 'bold';
            app.PlotMultipleGenotypesButton.Position = [133 101 187 32];
            app.PlotMultipleGenotypesButton.Text = 'Plot Multiple Genotypes';

            % Create ColumnsEditFieldLabel
            app.ColumnsEditFieldLabel = uilabel(app.PlottingTab);
            app.ColumnsEditFieldLabel.HorizontalAlignment = 'right';
            app.ColumnsEditFieldLabel.Position = [174 18 63 22];
            app.ColumnsEditFieldLabel.Text = '# Columns';

            % Create numColumns
            app.numColumns = uieditfield(app.PlottingTab, 'numeric');
            app.numColumns.Tooltip = {'Hint - set to zero to have the same number of columns as genotypes. '};
            app.numColumns.Position = [249 18 31 22];

            % Create PlotSettingsTab
            app.PlotSettingsTab = uitab(app.TabGroup);
            app.PlotSettingsTab.Title = 'Plot Settings';

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
            app.bulkYLim.Tooltip = {'Y Axis Limits [Min Max] (arbitrary units) for Bulk Signal plots.'};
            app.bulkYLim.Position = [98 8 74 22];
            app.bulkYLim.Value = '-500 5000';

            % Create ColorsPanel
            app.ColorsPanel = uipanel(app.PlotSettingsTab);
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

            % Create XAxesPanel
            app.XAxesPanel = uipanel(app.PlotSettingsTab);
            app.XAxesPanel.Title = 'X Axes';
            app.XAxesPanel.FontWeight = 'bold';
            app.XAxesPanel.Position = [24 239 174 84];

            % Create BulkSignallimitsLabel
            app.BulkSignallimitsLabel = uilabel(app.XAxesPanel);
            app.BulkSignallimitsLabel.HorizontalAlignment = 'right';
            app.BulkSignallimitsLabel.Position = [2 36 90 22];
            app.BulkSignallimitsLabel.Text = 'Bulk Signal lims';

            % Create bulkXLim
            app.bulkXLim = uieditfield(app.XAxesPanel, 'text');
            app.bulkXLim.Tooltip = {'X axis limits (minutes) for bulk signal traces.'};
            app.bulkXLim.Position = [130 36 39 22];
            app.bulkXLim.Value = '0 10';

            % Create BulkSignallimitsLabel_2
            app.BulkSignallimitsLabel_2 = uilabel(app.XAxesPanel);
            app.BulkSignallimitsLabel_2.HorizontalAlignment = 'right';
            app.BulkSignallimitsLabel_2.Position = [3 10 99 22];
            app.BulkSignallimitsLabel_2.Text = 'Axial Signal Ticks';

            % Create axialXTickInt
            app.axialXTickInt = uieditfield(app.XAxesPanel, 'text');
            app.axialXTickInt.Tooltip = {'Increments for X Axis Ticks (minutes) for Axial Signal plots'};
            app.axialXTickInt.Position = [130 10 39 22];
            app.axialXTickInt.Value = '2';

            % Create SortTracesByButtonGroup
            app.SortTracesByButtonGroup = uibuttongroup(app.PlotSettingsTab);
            app.SortTracesByButtonGroup.Title = 'Sort Traces By:';
            app.SortTracesByButtonGroup.FontWeight = 'bold';
            app.SortTracesByButtonGroup.Position = [216 194 120 103];

            % Create dontSort
            app.dontSort = uiradiobutton(app.SortTracesByButtonGroup);
            app.dontSort.Text = 'Dont sort';
            app.dontSort.Position = [13 50 71 22];
            app.dontSort.Value = true;

            % Create sortFreq
            app.sortFreq = uiradiobutton(app.SortTracesByButtonGroup);
            app.sortFreq.Text = 'Spike frequency';
            app.sortFreq.Position = [13 28 108 22];

            % Create sortAmp
            app.sortAmp = uiradiobutton(app.SortTracesByButtonGroup);
            app.sortAmp.Text = 'Spike amplitude';
            app.sortAmp.Position = [13 6 107 22];

            % Create HistogramBinsPanel
            app.HistogramBinsPanel = uipanel(app.PlotSettingsTab);
            app.HistogramBinsPanel.Title = 'Histogram Bins';
            app.HistogramBinsPanel.FontWeight = 'bold';
            app.HistogramBinsPanel.Position = [349 185 100 112];

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
            app.peakthreshold.Value = 500;

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
            app.FrameRateEditFieldLabel.Position = [347 67 68 22];
            app.FrameRateEditFieldLabel.Text = 'Frame Rate';

            % Create FrameRateEditField
            app.FrameRateEditField = uieditfield(app.PlotSettingsTab, 'numeric');
            app.FrameRateEditField.Tooltip = {'Frames per second of recorded data'};
            app.FrameRateEditField.Position = [420 69 29 19];
            app.FrameRateEditField.Value = 15;

            % Create NormalizeCheckBox
            app.NormalizeCheckBox = uicheckbox(app.PlotSettingsTab);
            app.NormalizeCheckBox.ValueChangedFcn = createCallbackFcn(app, @NormalizeCheckBoxValueChanged, true);
            app.NormalizeCheckBox.Text = 'Normalize Bulk Signal';
            app.NormalizeCheckBox.Position = [26 64 154 22];

            % Create AutoFixAxialSignalCheckBox
            app.AutoFixAxialSignalCheckBox = uicheckbox(app.PlotSettingsTab);
            app.AutoFixAxialSignalCheckBox.Tooltip = {'Re-analyze axial signal to correct for head tail flips.'};
            app.AutoFixAxialSignalCheckBox.Text = 'Auto-Fix Axial Signal';
            app.AutoFixAxialSignalCheckBox.Position = [26 28 132 24];
            app.AutoFixAxialSignalCheckBox.Value = true;

            % Create test
            app.test = uibutton(app.PlotSettingsTab, 'push');
            app.test.ButtonPushedFcn = createCallbackFcn(app, @testButtonPushed, true);
            app.test.Position = [229 64 82 43];

            % Create toQuerry
            app.toQuerry = uieditfield(app.PlotSettingsTab, 'numeric');
            app.toQuerry.Tooltip = {'fraction of axial signal to compare for head/tail reassignment'};
            app.toQuerry.Position = [24 7 42 19];
            app.toQuerry.Value = 0.15;

            % Create FractiontoQuerrylengthLabel
            app.FractiontoQuerrylengthLabel = uilabel(app.PlotSettingsTab);
            app.FractiontoQuerrylengthLabel.HorizontalAlignment = 'center';
            app.FractiontoQuerrylengthLabel.WordWrap = 'on';
            app.FractiontoQuerrylengthLabel.Position = [66 7 118 19];
            app.FractiontoQuerrylengthLabel.Text = '% of length to querry';

            % Create SortDirectionButtonGroup
            app.SortDirectionButtonGroup = uibuttongroup(app.PlotSettingsTab);
            app.SortDirectionButtonGroup.Title = 'Sort Direction';
            app.SortDirectionButtonGroup.FontWeight = 'bold';
            app.SortDirectionButtonGroup.Position = [216 126 120 62];

            % Create descendButton
            app.descendButton = uiradiobutton(app.SortDirectionButtonGroup);
            app.descendButton.Text = 'descend';
            app.descendButton.Position = [10 22 67 22];
            app.descendButton.Value = true;

            % Create ascendButton
            app.ascendButton = uiradiobutton(app.SortDirectionButtonGroup);
            app.ascendButton.Text = 'ascend';
            app.ascendButton.Position = [10 0 65 22];

            % Create BulkSignalPlotsButtonGroup
            app.BulkSignalPlotsButtonGroup = uibuttongroup(app.PlotSettingsTab);
            app.BulkSignalPlotsButtonGroup.Tooltip = {'should bulk signal be plotted in its own axis or overlaid on the axialSignal'};
            app.BulkSignalPlotsButtonGroup.Title = 'Bulk Signal Plots';
            app.BulkSignalPlotsButtonGroup.Position = [350 107 100 70];

            % Create SeparateButton
            app.SeparateButton = uiradiobutton(app.BulkSignalPlotsButtonGroup);
            app.SeparateButton.Text = 'Separate';
            app.SeparateButton.Position = [11 24 71 22];
            app.SeparateButton.Value = true;

            % Create OverlayButton
            app.OverlayButton = uiradiobutton(app.BulkSignalPlotsButtonGroup);
            app.OverlayButton.Text = 'Overlay';
            app.OverlayButton.Position = [11 2 65 22];

            % Create PlotsEditFieldLabel
            app.PlotsEditFieldLabel = uilabel(app.PlotSettingsTab);
            app.PlotsEditFieldLabel.HorizontalAlignment = 'right';
            app.PlotsEditFieldLabel.Position = [368 28 42 22];
            app.PlotsEditFieldLabel.Text = '# Plots';

            % Create numPlotsEditField
            app.numPlotsEditField = uieditfield(app.PlotSettingsTab, 'numeric');
            app.numPlotsEditField.Tooltip = {'number of bulk/axial signal plots to display. does not affect other graphs. '};
            app.numPlotsEditField.Position = [420 28 29 22];
            app.numPlotsEditField.Value = 10;

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