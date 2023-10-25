function save_struct_to_openPMD(gdfFilePath, openPMDOutputPath)
    % Define paths and names
    venvName = 'gdf_to_openpmd';
    venvPath = fullfile(pwd, venvName);
    pythonExe = fullfile(venvPath, 'bin', 'python'); % For macOS and Linux
    % pythonExe = fullfile(venvPath, 'Scripts', 'python.exe'); % Uncomment for Windows

    % Check if virtual environment already exists
    if ~isfolder(venvPath)
        % Create a new virtual environment
        system(sprintf('python -m venv %s', venvName));

        % Activate the virtual environment and install packages
        system(sprintf('%s -m pip install openpmd-api~=0.13.0,~=0.14.0, matplotlib>=3.0.0', pythonExe));
    end

    % Set MATLAB to use the Python from the virtual environment
    pyversion(pythonExe);

    % Call the gdf_to_openPMD.py script
    pythonScriptPath = 'gdf_to_openPMD.py'; % Update this to the actual path
    system(sprintf('%s %s -gdf %s -openPMD_output %s', pythonExe, pythonScriptPath, gdfFilePath, openPMDOutputPath));

    disp('Conversion completed!');
end
