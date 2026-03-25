d = dir('Z:\OAS\Patch_Foraging\**\*wormdata.mat');
int9pooled = [];
int1pooled = [];
for i =1:length(d)
    if i ~= 3
        fn = fullfile(d(i).folder, d(i).name)
        localFile = ['C:\tmp\' d(i).name];
        copyfile(fn, localFile);
        load(localFile)

        if isfield(wormdata, 'onFood')
            [int1, int9] = alignOffFood(wormdata);
            saveas(gcf, strrep(fn, 'wormdata.mat', 'exitTraces.png'))
            int1pooled  = horzcat(int1pooled, int1);
            int9pooled  = horzcat(int9pooled, int9);
            title(gca, i)
        end

        delete(localFile)
    end
end


