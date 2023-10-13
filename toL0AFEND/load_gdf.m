% Load binary gdf file from GPT in matlab structure array
%
% Input
%       filename    File to load
%       arrays_to_load      (Optional) arrays to load from the gdf file,
%                   specify as cell array {'name1';'name2'} etc
%                   For example {'x';'Bx'}
%                   All arrays are loaded when this option is not given,
%                   or when '' is passed. 
%       elements_to_load    (Optional) elements to load from the gdf file, 
%                   specified as an array. For example when you have a gdf
%                   file with particle output at t=0 and t=1 as a function
%                   of a scanned parameters, and you are only interested in
%                   t=1 for now, specify elements = [2 4 6 ...]
%                   When you are using this option and want to load all
%                   array elements (see above), pass '' as the second
%                   argument.
%                   When not given, all elements are loaded
%
% Output
%       gdf         1 x N structure array with data with N elements and 
%                   two fields for each element:
%                   .d  structure with data in single group. Fields:
%                     ID  Particle identification numbers
%                     x   x coordinate [m]
%                     y   y coordinate [m]
%                     z   z coordinate [m]
%                     Bx  Normalized velocity beta_x = v_x / c
%                     By  Normalized velocity beta_y = v_y / c
%                     Bz  Normalized velocity beta_z = v_y / c
%                     G   Lorentz factor
%                     rxy Distance to z-axis [m], sqrt(x^2+y^2)
%                     When using screen output command in GPT:
%                     t   Time when particles cross the screen [s] 
%                     When using tout (time output) command in GPT:
%                     fEx Ex field at the particle coordinates [V/m]
%                     fEy Ey field at the particle coordinates [V/m]
%                     fEz Ez field at the particle coordinates [V/m]
%                     fBx Bx field at the particle coordinates [T]
%                     fBy By field at the particle coordinates [T]
%                     fBz Bz field at the particle coordinates [T]
%                   .p  structure with params of that group
%       info        structured array with main header information
%
% Example
%       load_gdf('data.gdf')           %Load whole file
%       load_gdf('data.gdf',{'x';'y'}) %Load only x and y data
%       load_gdf('data.gdf','', [2 4]) %Load only elements 2 and 4
%
% Structure gdf file: (each line contains a block of data)
% Header (level = 1 after this step)
% param @logo
% params scanned by mr (level = 2 after this step)
% param creator
% param @logo
% param position / time (level = 3 after this step)
% dat particle data (level = 2 after this step)
% param position / time (possibly) (level = 3 after this step)
% dat particle data     (possibly) (level = 2 after this step)
% param numderivs (only loaded into structure when elements_to_load is not specified as input)
% param cputime   (only loaded into structure when elements_to_load is not specified as input)
% (level = 1 after this step). All data belonging to these mr settings are read
% params scanned by mr
% param @logo
% param position / time
% dat particle data
% param position / time (possibly)
% dat particle data     (possibly)
% param numderivs (only loaded into structure when elements_to_load is not specified as input)
% param cputime   (only loaded into structure when elements_to_load is not specified as input)
% etc
%
%
function [gdf, info] = load_gdf(filename, arrays_to_load, elements_to_load)

if nargin == 1 
    arrays_to_load   = []; %load all arrays
    elements_to_load = []; %load all elements
    last_element_to_load = -1; %load all elements
elseif nargin == 2 
    elements_to_load = []; %load all elements
    last_element_to_load = -1; %load all elements
else
    last_element_to_load = max(elements_to_load); % stop when this element is loaded
end

% Constants
GDFNAMELEN = 16;                 % Length of the ascii-names
GDFID  = 94325877;               

% Data types 
t_ascii  = hex2dec('0001');      % Char array
t_s32    = hex2dec('0002');      % Signed long           
t_dbl    = hex2dec('0003');      % Double    
t_nul    = hex2dec('0010');      % No data               

% t_undef  = hex2dec('0000');      % Data type not defined 
% t_u8     = hex2dec('0020');      % Unsigned char         
% t_s8     = hex2dec('0030');      % Signed char           
% t_u16    = hex2dec('0040');      % Unsigned short        
% t_s16    = hex2dec('0050');      % Signed short          
% t_u32    = hex2dec('0060');      % Unsigned long         
% t_u64    = hex2dec('0070');      % Unsigned 64bit int    
% t_s64    = hex2dec('0080');      % Signed 64bit int      
% t_flt    = hex2dec('0090');      % Float                 

% Block types 
t_dir    = hex2dec('0100');      % Directory entry start
t_edir   = hex2dec('0200');      % Directory entry end 
t_param  = hex2dec('0400');      % Parameter        
t_data   = hex2dec('0800');      % Data array

% Open file
FileID = fopen(filename, 'r', 'l', 'windows-1252'); % Specify the encoding of the file, otherwise this function won't work under linux
if (FileID<0)
    error('Error in ''load_gdf'': did not find file ''%s''', filename);
end

% Read gdfmainheader
ID = fread(FileID, 1, '*uint32');
if (ID ~= GDFID) % check if it is indeed a GDF file
    fclose(FileID); % close file
	error('Error in ''load_gdf'': this is not a gdf file');    
end
info = struct;
info.ID = ID;
info.cretime  = fread(FileID, 1, '*uint32');
info.creator  = fread(FileID, GDFNAMELEN, '*char')';
info.creator  = info.creator(1: find( (info.creator==0),1)-1);      % get string part upto zero-character
info.destin   = fread(FileID, GDFNAMELEN, '*char')';
info.destin   = info.destin(1: find( (info.destin==0),1)-1);        % get string part upto zero-character
info.gdfmaj   = fread(FileID, 1, '*uint8');
info.gdfmin   = fread(FileID, 1, '*uint8');
info.cremaj   = fread(FileID, 1, '*uint8');
info.cremin   = fread(FileID, 1, '*uint8');
info.desmaj   = fread(FileID, 1, '*uint8');
info.desmin   = fread(FileID, 1, '*uint8');
dummy         = fread(FileID, 2, '*uint8');


% Read gdf data blocks
gdf    = struct; % structure which will be filled with data from GDF
params = struct; % contains the scanned parameters when using MR, time when using tout, position for screen
params.level1 = [];
arrays = struct; % contains particle info: position, momenta, electromagnetic fields
current_element = 1;
elements_stored = 0; % contains the number of elements that are already saved
level = 1;

while (~feof(FileID))    
    %read gdf block header
    name    = fread(FileID, GDFNAMELEN, '*char')';
    name    = name(1: find( (name==0),1)-1);      % get string part up to zero-character
    
    % get rid of not allowed characters etc
    if isempty(name)
        name = 'empty';
    elseif name(1) == '@'
       name = 'x0x40logo';        
    end
    
    block_type = fread(FileID, 1, '*uint32');        % block type
    data_type  = bitand(block_type, 255);            % get byte with datatype
    block_size = fread(FileID, 1, '*uint32');        % byte count of data in block
    
    % get block type
    start_dir     = ( bitand(block_type, t_dir)  >0);
    end_dir       = ( bitand(block_type, t_edir) >0); % end folder 
    data_is_param = ( bitand(block_type, t_param)>0); % if 1, the data that is coming next in the file contains the scanned parameters when using MR, time when using tout, position for screen
    data_is_array = ( bitand(block_type, t_data) >0); % if 1, the data that is coming next in the file contains particle info: position, momenta, electromagnetic fields
    
    % for debugging
%     fprintf('%-16s', name)
%     if start_dir
%         fprintf('\t# Start dir. Level = %d', level+1)
%     end
%     if end_dir
%         fprintf('\t# End dir.   Level = %d', level-1)
%     end
%     fprintf('\n')
    
    % new folder
    if (start_dir)
        level = level + 1;
        % Parameters for the new level are parameters for one level lower
        % plus parameters that will be loaded for this level
        % Example: when scanning with MR, the parameters are defined in
        % level 2, while the time/position output is defined in level 3. We
        % need to copy this info, so both the scanned parameter and the
        % time/output info is loaded for this data group
        params.(sprintf('level%d',level)) = params.(sprintf('level%d',level-1));
    end
    
    % end folder 
    if end_dir % when end_dir = 1, data_is_param = 1, name = empty, and data_type = t_nul (no data)                                         
        if isempty(fieldnames(arrays))
            % params contains numderivs and cpu. 
            if isempty(elements_to_load)
                % Attach to last data point  
                if ~iscell(arrays_to_load) && strcmp(arrays_to_load, ' ')
                    % no array data is loaded. Store only parameter data
                    elements_stored = elements_stored + 1;
                elseif  elements_stored == 0
                     elements_stored = 1;
                end
                gdf(elements_stored).p = params.(sprintf('level%d',level));
            end
            % To implement: save into element when not all elements are
            % loaded. Challenge: you have to figure out if this data
            % belongs to one of the saved elements
        else
            % Check if this element should be saved
            if (isempty(elements_to_load)) || any(elements_to_load == current_element)
                % save data group
                elements_stored = elements_stored + 1;
                gdf(elements_stored).d = arrays;
                gdf(elements_stored).p = params.(sprintf('level%d',level));
                if current_element == last_element_to_load
                    % last requested element loaded. Close and exit
                    break        
                end
            end
            arrays = struct;                         % clear data arrays
            current_element = current_element + 1;   % next data group  
        end
        level = level - 1; 
    end   

    % read a scanned parameters when using MR, time when using tout, position for screen
    if data_is_param      
        switch data_type
            case t_dbl
                value = fread(FileID, 1, 'double');                   %read param
                params.(sprintf('level%d',level)).(name) = value;     %save param
            case t_nul
                % no data present
            case t_ascii
                value  = fread(FileID, double(block_size), '*char')'; %read param
                params.(sprintf('level%d',level)).(name) = value;     %save param
            case t_s32
                value  = fread(FileID, 1, 'int32')';                  %read param
                params.(sprintf('level%d',level)).(name) = value;     %save param                
            otherwise
                % error, abort while loop
                fclose(FileID); % close file
                error('Error in ''load_gdf'': unknown datatype of value\nByte position in file: %d\n', ftell(FileID))
        end
    end    
    
    % read data array with particle info: position, momenta, electromagnetic fields
    if data_is_array
        switch data_type
            case t_dbl
                if (mod(block_size,8) ~= 0); 
                    fclose(FileID); % close file
                    error('Error in ''load_gdf'': wrong size for double array'); 
                end;
                N = double(block_size/8);                             %number of array elements
                value = fread(FileID, N, 'double');                   %read array

                if ( isempty(arrays_to_load) ) || any( strcmp(name, arrays_to_load) ) % check if array should be saved
                    try
                        arrays.(name) = value; 
                    catch
                        % There are some non allowed characters in name. Get rid of them
                        name = genvarname(name);
                        arrays.(name) = value;
                    end    
                end;
            otherwise
                % error, abort while loop
                fclose(FileID); % close file
                error('Error in ''load_gdf'': unknown datatype of array\nByte position in file: %d\n', ftell(FileID))
        end
        
    end       
end

% params contains numderivs and cpu for the last element. 
% these are not added to the structure yet. Do this here
if isempty(elements_to_load) && isfield(params, sprintf('level%d',level) ) && ...
        isfield(params.(sprintf('level%d',level)), 'cputime');
    % Attach to last data point  
    if  elements_stored == 0
         elements_stored = 1;
    end
    level_str = sprintf('level%d',level);
    field_names = fieldnames( params.(level_str) );
    for i=1:length(field_names)
        gdf(elements_stored).p.(field_names{i}) = params.(level_str).(field_names{i});
    end
end
% To implement: save into element when not all elements are
% loaded. Challenge: you have to figure out if this data
% belongs to one of the saved elements


% Some gdf files do not use start / end dir identifiers. Save the data here
% then
if ~isempty(fieldnames(arrays)) && isempty(elements_to_load)
    % save data group
    elements_stored = elements_stored + 1;
    gdf(elements_stored).d = arrays;
    gdf(elements_stored).p = params.(sprintf('level%d',level));
end

fclose(FileID); % close file

end

