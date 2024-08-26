function [] = fixAxialSignal()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
[fn,fp] = uigetfile("Z:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\Receptor_Mutants\wildtype-30mM-TA\220217_zfis178_wildtype-30mM-TA_1");
path = fullfile(fp,fn);
tempdir = 'C:\tmp';
tempfullfile = fullfile(tempdir, fn);


if isfile(tempfullfile) == 0
    disp(['Copying File: ' tempfullfile ' to temp dir'])
    copyfile(path, tempdir);
end

m = matfile(tempfullfile, 'Writable',true);
data = m.wormdata;
data.temppath = tempfullfile;
data.serverpath = path;

% if ~isfield(data, 'fixedSignal')
%     data.workingSignal = data.rawAxialSignal;
% elseif isfield(data, 'fixedSignal')
%     data.workingSignal = data.fixedSignal;
%     disp('Using previously modified signal');
% end
data.workingSignal = data.autoAxialSignal;

% data.outname = axsigname;
fig = figure('Position', [506.6000 33 560.0000 749.6000]);
uax = axes('Parent',fig,'Position', ...
    [0.1300 0.03100 0.50 0.950]);
imagesc(smoothdata(data.workingSignal,'gaussian',60))

colormap(uax, viridis)
btn = uicontrol('Style', 'pushbutton',...
    'Position',[400, 700, 120, 22],'String', 'Mark Points',...
    'Callback', @btn_callback);

undo = uicontrol('Style', 'pushbutton',...
    'Position',[400, 650, 120, 22],'String', 'Undo/Redo',...
    'Callback', @undo_callback);

flpallbtn = uicontrol('Style', 'pushbutton',...
    'Position',[400, 625, 120, 22],'String', 'Flip Entire Signal',...
    'Callback', @flpallbtn_callback);

rvert = uicontrol('Style', 'pushbutton',...
    'Position',[400, 575, 120, 22],'String', 'Revert to Raw Signal',...
    'Callback', @rvert_callback);


loadfixed = uicontrol('Style', 'pushbutton',...
    'Position',[400, 550, 120, 22],'String', 'Revert to Fixed Signal',...
    'Callback', @loadfixed_callback);


autofix = uicontrol('Style', 'pushbutton',...
    'Position',[400, 500, 120, 22],'String', 'Attempt Auto Fix',...
    'Callback', @autofix_callback);

endbtn = uicontrol('Style', 'pushbutton',...
    'Position',[400, 450, 120, 22],'String', 'Save Signal',...
    'Callback', @endbtn_callback);

cleanup = uicontrol('Style', 'pushbutton',...
    'Position',[400, 425, 120, 22],'String', 'clear and upload',...
    'Callback', @cleanup_callback);

guidata(fig, data);


    function btn_callback(fig, ~)
        fig = gcf;
        data= guidata(fig);
        [~,ypts]= ginput(2);
        ystart = floor(ypts(1));
        yend = floor(ypts(2));
        
        data.workingSignal(ystart:yend,:) = fliplr(data.workingSignal(ystart:yend,:));
        imagesc(smoothdata(data.workingSignal,'gaussian',60)) 
        
        colormap(uax, viridis)
        
        if ~isfield(data, 'ypts')
            data.ypts = [ystart yend];
        else
            data.ypts = vertcat(data.ypts, ypts');
        end
        
        data.ypts
        
        guidata(fig, data);
        drawnow
    end

    function undo_callback(fig, ~)
        data= guidata(fig);
        ystart = floor(data.ypts(end,1));
        yend = floor(data.ypts(end,2));
        data.workingSignal(ystart:yend,:) = fliplr(data.workingSignal(ystart:yend,:));
        imagesc(smoothdata(data.workingSignal,'gaussian',60)) 
        
        colormap(uax, viridis)
        guidata(fig, data);
    end

    function flpallbtn_callback(fig, ~)
        data= guidata(fig);
        ypts = [1 length(data.workingSignal)];
        data.workingSignal(ypts(1):ypts(2),:) = fliplr(data.workingSignal(ypts(1):ypts(2),:));
        imagesc(smoothdata(data.workingSignal,'gaussian',60)) 
        
        colormap(uax, viridis)
        
        data.ypts = ypts;
        guidata(fig, data);
        drawnow
    end

    function rvert_callback(fig, ~)
        data= guidata(fig);
        data.workingSignal = data.rawAxialSignal;
        imagesc(smoothdata(data.workingSignal,'gaussian',60)) 
        
        colormap(uax, viridis)
        
        guidata(fig, data);
    end

    function loadfixed_callback(fig, ~)
        data= guidata(fig);
        
        if isfield(data, 'fixedSignal')
            data.workingSignal = data.fixedSignal;
        else
            disp('Sorry, no fixed signal found');
        end
        
        imagesc(smoothdata(data.workingSignal,'gaussian',60)) 
        
        colormap(uax, viridis)
        
        guidata(fig, data);
    end


    function autofix_callback(fig,~)
        data = guidata(fig);
        for i = 1:length(data.workingSignal)
            left = mean(data.workingSignal(i,1:20),'omitnan');
            right = mean(data.workingSignal(i,end-20:end),'omitnan');
            if left > right
                data.workingSignal(i,:) = fliplr(data.workingSignal(i,:));
            end
        end
        imagesc(smoothdata(data.workingSignal,'gaussian',60)) 
        colormap(uax, viridis)
        
        guidata(fig, data);
    end


    function endbtn_callback(fig,~)
        data = guidata(fig);
        data.fixedSignal = data.workingSignal;
        load(data.temppath)
        wormdata.autoAxialSignal = data.workingSignal;
        wormdata.noAutoFix = 1;
        
        if isfield(data, 'ypts')
            wormdata.axialTransform = data.ypts;
        end


        save(data.temppath, 'wormdata')
        guidata(fig, data);
        
        disp('data saved locally, ready to upload');
    end

    function cleanup_callback(fig,~)
        data = guidata(fig);
        disp(['moving local data to: ' data.serverpath])
        try
            [status, message, ID] = copyfile(data.temppath, data.serverpath);
            if status == 0
                disp(message)
            elseif status == 1
                
                disp('deleting local data');
                delete(data.temppath);
            end
        catch
        end
    end




end

