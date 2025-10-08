function params = loadHHGPMParameters(parameterFilename)
%LOADHHGPM_PARAMETERS Read the text parameter file used by HHG_Gouy_PM_v4.
%   The legacy implementation relied on executing each line via EVAL. The
%   refactored version keeps the same flexibility but stores the results in
%   a structure so that downstream routines can explicitly request the
%   values they need.
%
%   PARAMS = LOADHHGPM_PARAMETERS(FILENAME) returns a structure whose field
%   names match the variable names defined inside the parameter file.

   fid = fopen(parameterFilename, 'r');
   if fid < 0
       error('loadHHGPMParameters:FileNotFound', ...
             'Unable to open parameter file "%s".', parameterFilename);
   end

   cleaner = onCleanup(@() fclose(fid));

   params = struct();
   while ~feof(fid)
       line = fgetl(fid);
       if ~ischar(line)
           break;
       end

       line = strtrim(line);
       if isempty(line) || startsWith(line, '%')
           continue;
       end

       % Prepend the structure name so that assignments are stored inside
       % the params structure. Lines in the parameter file end with a
       % semicolon, therefore no extra punctuation is required here.
       eval(sprintf('params.%s', line)); %#ok<EVLP>
   end
end
