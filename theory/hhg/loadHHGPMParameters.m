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

   % Define the unit-conversion helpers that existed in the legacy script so
   % that parameter files which reference them (for example ``808*nm``) keep
   % working when parsed in this sandboxed scope.
   eV = 1.602e-19; %#ok<NASGU>
   cm = 1e-2;      %#ok<NASGU>
   mm = 1e-3;      %#ok<NASGU>
   um = 1e-6;      %#ok<NASGU>
   nm = 1e-9;      %#ok<NASGU>
   fs = 1e-15;     %#ok<NASGU>
   mJ = 1e-3;      %#ok<NASGU>

   c = 2.99792458e8; %#ok<NASGU>  % speed of light (m/s)
   q_e = 1.602e-19;   %#ok<NASGU>  % electron charge (C)
   m_e = 9.101e-31;   %#ok<NASGU>  % electron mass (kg)
   mu_0 = 4*pi*1e-7;  %#ok<NASGU>  % permeability of vacuum (N/A^2)
   epsilon_0 = 8.854e-12; %#ok<NASGU> % permittivity of vacuum (F/m)
   h = 6.626e-34;     %#ok<NASGU>  % Planck constant (J*s)
   hbar = h/(2*pi);   %#ok<NASGU>  % Reduced Planck constant (J*s)

   
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
