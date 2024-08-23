% SAVEBADNAMES saves a mat file, using intentionally bad names as it goes
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    savebadnames
%  Filename:    savebadnames.m
%  Programmer:  James Tursa
%  Version:     1.00
%  Date:        June 18, 2013
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
%  SAVEBADNAMES is typically self building. That is, the first time you call it,
%  the savebadnames.m file recognizes that the mex routine needs to be compiled and
%  then the compilation will happen automatically.
% 
%  The usage is as follows (arguments in brackets [ ] are optional):
% 
%  Syntax
% 
%  savebadnames(FILENAME [,verbose])
% 
%      FILENAME = The mat file to be saved
%      verbose  = 1 or 0 (optional, causes name change log to be displayed)
% 
%  Description
% 
%  SAVEBADNAMES saves a mat file using intentionally invalid names.
%  This is a program to create a test mat file for the LOADFIXNAMES function.
% 
%  Change Log:
%  2013/Jun/18 --> 1.00, Initial Release
%--
%
%**************************************************************************

function varargout = savebadnames(varargin)
disp(' ');
disp('... Mex routine has not been built yet ... building it ...');

%\
% Check to see that loadfixnames.c source code is present
%/

disp('... Finding path of savebadnames C source code files');
try
    mname = mfilename('fullpath');
catch
    mname = mfilename;
end
cname = [mname '.c'];
if( isempty(dir(cname)) )
    disp('Cannot find the file savebadnames.c in the same directory as the');
    disp('file savebadnames.m. Please ensure that they are in the same');
    disp('directory and try again. The following file was not found:');
    disp(' ');
    disp(cname);
    disp(' ');
    error('Unable to compile savebadnames.c');
end
disp(['... Found file savebadnames.c in ' cname]);

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
    disp('... mex savebadnames.c build completed ... you may now use savebadnames.');
    disp(' ');
catch
    disp(' ');
    disp('... Well, *that* didn''t work!');
    disp('... Please contact author ...');
    disp(' ');
    disp('... Changing to original directory ...');
    cd(cdold);
    error('Unable to compile savebadnames.c');
end

%\
% Restore old directory
%/

disp('... Changing to original directory ...');
cd(cdold);

%\
% Call the mex routine savebadnames.
%/

[varargout{1:nargout}] = savebadnames(varargin{:});

return
end
