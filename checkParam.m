function [param] = checkParam(param)
% Check if the parameters are assigned reasonable values. If no value is
% assigned by the user, the default value will be used for that parameter.
% Author: Rizwan Ahmad (ahmad.46@osu.edu)


% param.typ
if ~isfield(param, 'typ') || isempty(param.typ)
    error('No value has been assigned to param.typ');
elseif ~strcmp(param.typ, 'UIS') && ~strcmp(param.typ, 'VRS') && ~strcmp(param.typ, 'VISTA')
    error('param.typ can only assume one of the three values: UIS, VRS, or VISTA');
end

% param.p
if ~isfield(param, 'p') || isempty(param.p)
    error('No value has been assigned to param.p');
elseif param.p < 2 || rem(param.p,1) ~= 0
    error('The value assigned to param.p must be an integer greater than 1');
end

% param.t
if ~isfield(param, 't') || isempty(param.t)
    error('No value has been assigned to param.t');
elseif param.t < 2 || rem(param.t,1) ~= 0
    error('The value assigned to param.t must be an integer greater than 1');
elseif param.t > 32 && strcmp(param.typ, 'VISTA');
    fprintf(['----------------------------------- Tip ------------------------------------' '\n' ...
             'For faster processing, reduce the number of frames (to lets say 32) and then' '\n' ...
             'cyclically reuse these frames to achieve any arbitrarilty number of frames.' '\n'...
             '----------------------------------------------------------------------------' '\n']);
end

% param.R
if ~isfield(param, 'R') || isempty(param.R)
    error('No value has been assigned to param.R');
elseif param.R < 1 || param.R > param.p
    error('The value assigned to param.R must be an integer between 1 and param.p');
end

% param.alph
if ~isfield(param, 'alph') || isempty(param.alph)
    error('No value has been assigned to param.alph');
elseif param.alph < 0 || param.alph > 1
    error('The value assigned to param.alph must be between 0 and 1');
end

% param.sig
if ~isfield(param, 'sig') || isempty(param.sig)
    error('No value has been assigned to param.sig');
elseif param.sig <= 0
    error('The value assigned to param.sig must be greater than zero');
end


%% Parameters with preset default values
% param.nIter
if ~isfield(param, 'nIter') || isempty(param.nIter)
    param.nIter = 120;
elseif param.nIter < 20 || param.nIter > 1024 || rem(param.nIter,1) ~= 0
    error('The value assigned to param.nIter must be an integer between 20 and 1024');
end

% param.ss
if ~isfield(param, 'ss') || isempty(param.ss)
    param.ss = 0.25; % default
elseif param.ss <= 0
    error('The value assigned to param.ss must be greater than zero');
end

% param.tf
if ~isfield(param, 'tf') || isempty(param.tf)
    param.tf = 0.0; % default
elseif param.tf < 0
    error('The value assigned to param.tf must be greater than or equql to zero');
elseif param.tf > 0
    warning('param.tf > 0 will destroy the constant temporal resolution.');
end

% param.s
if ~isfield(param, 's') || isempty(param.s)
    param.s = 1.4; % default
elseif param.s <= 0 || param.s > 10
    error('The value assigned to param.s must be between 0 and 10');
end

% param.g
if ~isfield(param, 'g') || isempty(param.g)
    param.g = floor(param.nIter/6); % default
elseif param.g < 5 || param.g > round(param.nIter/4) || rem(param.g,1) ~= 0
    error('The value assigned to param.g must be between 1 and param.nIter/4');
end

% param.uni
if ~isfield(param, 'uni') || isempty(param.uni)
    param.uni = round(param.nIter/2); % default
elseif param.uni < round(param.nIter/4) || param.uni > round(param.nIter/2) || rem(param.uni,1) ~= 0
    error('The value assigned to param.uni must be between param.nIter/4 and param.nIter/2');
end

% param.fl
if ~isfield(param, 'fl') || isempty(param.fl)
    param.fl= round(5/6*param.nIter); % default
elseif param.fl <= param.nIter/2 || param.fl > param.nIter || rem(param.fl,1) ~= 0
    error('The value assigned to param.fl must be between param.nIter/2 and param.nIter');
end

% param.W
if ~isfield(param, 'W') || isempty(param.W)
    param.W = max(param.R/10 + 0.25, 1);
%     param.W = param.R/10 + 0.25;
elseif param.W <= 0
    error('The value assigned to param.W must be greater than zero');
end

% param.sz
if ~isfield(param, 'sz') || isempty(param.sz)
    param.sz = 3.5; % default
elseif param.sz <= 0
    error('The value assigned to param.sz must be greater than zero');
end

% param.dsp
if ~isfield(param, 'dsp') || isempty(param.dsp)
    param.dsp = 1; % default
elseif param.dsp < 1 || param.dsp > param.nIter
    error('The value assigned to param.dsp must be between 1 and nIter');
end

% param.fs
if ~isfield(param, 'fs') || isempty(param.fs)
    param.fs = 1; % default
elseif (param.fs ~= 0 && param.fs ~= 1 && param.fs < param.R) || rem(param.fs,1)
    error('The value assigned to param.fs must be a non-negative integer');
end
