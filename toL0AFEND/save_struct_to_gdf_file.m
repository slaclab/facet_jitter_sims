% Save a matlab structure array to a binary gdf file for use in GPT
%
% Input
%       filename    File to save
%       data        Data structure to save
%
% Output
%       -
%
% Example
% data = struct;
% data.p.testdouble = 11;
% data.p.testascii = 'Test Ascii string';
% data.d.testx = 1:10;
% data.d.testy = 1:10;
% save_struct_to_gdf_file('test.gdf', data)
%
function save_struct_to_gdf_file(filename, data)

creator_string = 'MatlabGDFSave';

% Open file for writing
FileID = fopen(filename, 'w');
if FileID < 0
    error('Error: did not find file ''%s''', filename);
end

% Write GDF header
writegdf_mainhead(FileID, creator_string)

% Write data
for i=1:length(data)
    % Write parameters
    if isfield(data(i), 'p') && ~isempty(data(i).p)
        names = fieldnames(data(i).p);
        for j=1:length(names)
            writegdf_parameter(FileID,names{j},data(i).p.(names{j}))
        end
    end

    % Write arrays
    if isfield(data(i), 'd') && ~isempty(data(i).d)
        names = fieldnames(data(i).d);
        for j=1:length(names)
            writegdf_doublearray(FileID,names{j},data(i).d.(names{j}))
        end
    end

    % End this data block
    writegdf_endgroup(FileID);
end

fclose(FileID);

end


% Ensure GPT compatible name
function sname = writegdf_getname(s)

% Constants
GDFNAMELEN = 16;                 % Length of the ascii-names

if length(s) > GDFNAMELEN % Trim name, so maximum length is 16
    sname = s(1:16);
else
    sname = horzcat(s,zeros(1,GDFNAMELEN-length(s)));
end

end


% Write binary gdf file header
function writegdf_mainhead(FileID, creator)

% Constants
GDFNAMELEN = 16;                 % Length of the ascii-names
GDFID  = 94325877;               

fwrite(FileID, GDFID, '*uint32'); % GdfID
% The posixtime function is not available in Matlab 2012
% fwrite(FileID, posixtime(datetime('now')), '*uint32'); % Cretime
% Convert time to POSIX format
fwrite(FileID, (now-datenum('1970-1-1 00:00:00'))*(24*3600.0), '*uint32'); % Cretime
fwrite(FileID, writegdf_getname(creator), '*char'); % Creator
fwrite(FileID, zeros(1,GDFNAMELEN), '*char'); % Destin
fwrite(FileID, 1, '*uint8'); % Gdfmaj
fwrite(FileID, 1, '*uint8'); % Gdfmin
fwrite(FileID, 0, '*uint8'); % Cremaj
fwrite(FileID, 0, '*uint8'); % Cremin
fwrite(FileID, 0, '*uint8'); % Desmaj
fwrite(FileID, 0, '*uint8'); % Desmin
fwrite(FileID, 0, '*uint8'); % Dummy1
fwrite(FileID, 0, '*uint8'); % Dummy2

end


% Write an array to the gdf file
function writegdf_doublearray(FileID,name,data)

% Data types           
t_dbl    = hex2dec('0003');      % Double    
t_arr    = hex2dec('0800');      % Array  

% Write header
fwrite(FileID, writegdf_getname(name), '*char'); % name
fwrite(FileID, t_dbl + t_arr,'*uint32'); % block type
fwrite(FileID, 8*length(data),'*uint32'); % block size

% Write data
fwrite(FileID, data,'double');

end


% Write a parameter to the gdf file
function writegdf_parameter(FileID,name,data)

% Data types           
t_ascii  = hex2dec('0001');      % Char array
t_dbl    = hex2dec('0003');      % Double    
% t_s32    = hex2dec('0002');      % Signed long           
% signed long not yet implemented. For now use double only. There are a limited number of
% params writted as long (numderivs)
t_dir    = hex2dec('0100');      % Directory entry start
t_param  = hex2dec('0400');      % Parameter 

if isnumeric(data) % Data = double
    % Write header
    fwrite(FileID, writegdf_getname(name), '*char'); % name
    if strcmp(name, 'position') || strcmp(name, 'time') % start new data group
        fwrite(FileID, t_dbl + t_param + t_dir,'*uint32'); % block type
    else
        fwrite(FileID, t_dbl + t_param,'*uint32'); % block type
    end
    fwrite(FileID, 8,'*uint32'); % block size

    % Write data
    fwrite(FileID, data,'double');
else % Data = char array
    % Write header
    fwrite(FileID, writegdf_getname(name), '*char'); % name
    fwrite(FileID, t_ascii + t_param,'*uint32'); % block type
    fwrite(FileID, 1*length(data),'*uint32'); % block size

    % Write data
    fwrite(FileID, data,'*char');
end

end


% End a data block
function writegdf_endgroup(FileID)

name = '';

% Data types           
t_nul    = hex2dec('0010');      % No data
t_edir   = hex2dec('0200');      % Directory entry end 

% Write header
fwrite(FileID, writegdf_getname(name), '*char'); % name
fwrite(FileID, t_nul + t_edir,'*uint32'); % block type
fwrite(FileID, 0,'*uint32'); % block size

% No data to write

end
