% LOADFIXNAMES loads a mat file, fixing names as it goes
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    loadfixnames
%  Filename:    loadfixnames.m
%  Programmer:  James Tursa
%  Version:     1.20
%  Date:        October 31, 2013
%  Copyright:   (c) 2013 by James Tursa, All Rights Reserved
%
%  This code uses the BSD License:
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
%
%  Building:
% 
%  LOADFIXNAMES is typically self building. That is, the first time you call it,
%  the loadfixnames.m file recognizes that the mex routine needs to be compiled and
%  then the compilation will happen automatically.
% 
%  The usage is as follows (arguments in brackets [ ] are optional):
% 
%  Syntax
% 
%  loadfixnames(FILENAME [,names] [,verbose])
%  S = loadfixnames(FILENAME [,names] [,verbose])
% 
%      FILENAME = The mat file to be loaded
%      names    = string or cell array of strings (variable name(s) to load)
%      verbose  = 1 (optional, causes name change log to be displayed)
%      S = struct returned instead of loading variables into workspace
% 
%  Description
% 
%  LOADFIXNAMES loads a mat file into the workspace, fixing invalid names.
%  All invalid characters are replaced with an underscore. Also, if the
%  first character is not a letter, an 'A' is added to the front. If the
%  variable is a structure, fixes the field names also. In the case of
%  field names, the name length is kept constant. So if a field name
%  begins with a digit it will be replaced with 'A' - 'J' instead. If the
%  field name begins with an invalid non-digit it will be replaced with
%  'A' - 'Z' or 'a' - 'z' (letters are cycled in an attempt to avoid
%  name clashes).
% 
%  Limitations: The current renaming scheme makes a mild attempt to avoid
%               name clashes, but does not guarantee this.
%
%  See the companion function SAVEBADNAMES for a test function that will
%  intentionally create a mat file with invalid names.
% 
%  Change Log:
%  2013/Jun/06 --> 1.00, Initial Release to Answers (used deep copy)
%  2013/Jun/18 --> 1.10, Initial Release to FEX (uses shared data copy)
%  2013/Oct/31 --> 1.20, Added code to deal with duplicate variable names
%                        Added capability to return a struct
%--
%
%**************************************************************************

function varargout = loadfixnames(varargin)
disp(' ');
disp('... Mex routine has not been built yet ... building it ...');

%\
% Check to see that loadfixnames.c source code is present
%/

disp('... Finding path of loadfixnames C source code files');
try
    mname = mfilename('fullpath');
catch
    mname = mfilename;
end
cname = [mname '.c'];
if( isempty(dir(cname)) )
    disp('Cannot find the file loadfixnames.c in the same directory as the');
    disp('file loadfixnames.m. Please ensure that they are in the same');
    disp('directory and try again. The following file was not found:');
    disp(' ');
    disp(cname);
    disp(' ');
    error('Unable to compile loadfixnames.c');
end
disp(['... Found file loadfixnames.c in ' cname]);

%\
% Save old directory and change to source code directory
%/

disp('... Changing to source code directory ...');
cdold = cd;
if( length(mname) > 13 )
    cd(mname(1:end-13));
end

%\
% Do the compile
%/

disp('... Now attempting to compile ...');
disp(['mex(''' cname ''')']);
disp(' ');
try
    mex(cname);
    disp('... mex loadfixnames.c build completed ... you may now use loadfixnames.');
    disp(' ');
catch
    disp(' ');
    disp('... Well, *that* didn''t work!');
    disp('... Please contact author ...');
    disp(' ');
    disp('... Changing to original directory ...');
    cd(cdold);
    error('Unable to compile loadfixnames.c');
end

%\
% Restore old directory
%/

disp('... Changing to original directory ...');
cd(cdold);

%\
% Call the mex routine loadfixnames.
%/

disp(' ');
disp('Use up arrow to reissue the loadfixnames command.');
disp(' ');
if( 0 )
    [varargout{1:nargout}] = loadfixnames(varargin{:});
end

return
end
