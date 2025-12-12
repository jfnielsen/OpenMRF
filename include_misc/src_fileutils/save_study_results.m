function save_study_results(study_info, res, default_dir)
    savename = sprintf('%s.mat', study_info.study_name(1:end-4));
    if nargin < 3
        [filename, pathname] = uiputfile('*.mat', 'Save results as', savename);
    else 
        [filename, pathname] = uiputfile('*.mat', 'Save results as', savename, 'Folder', default_dir);
    end
    if nargin==1
        fprintf('no variables to save')
        return 
    end

    % If user didn't cancel, save to the chosen location
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('Save canceled by user.');
    else
        fullpath = fullfile(pathname, filename);
        save(fullpath, 'res');
        fprintf('Saved variables to: %s\n', fullpath);
    end
end