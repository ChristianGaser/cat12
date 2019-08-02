function cat_io_send_to_server
% This function will send Matlab version information to the SBM server
% The mac address is only used for creating a uniqe id to not count visits
% several times
% If you don't want to send this information to the SBM server 
% (for internal use only) change the flag cat.extopts.send_info in 
% cat_defaults.m
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

% don't do anything if default is not set
if ~cat_get_defaults('extopts.send_info'), return; end

mac = MACAddress();

% remove "-" or ":" in mac address
mac = strrep(mac,'-','');
mac = strrep(mac,':','');

% use 16 hex digits
mac = ['0000' mac];

url = sprintf('http://www.neuro.uni-jena.de/piwik/piwik.php?idsite=1&rec=1&_id=%s&action_name=%s%s%s',mac,'Start','%2F',version('-release'));
urlread(url);

function [mac, st] = MACAddress(allMac)
% [mac, st] = MACAddress()
% 
% The default is to return one MAC address, likely for ethernet adaptor. If the
% optional input is provided and true, all MAC address are returned in cellstr.
% No internet connection is required for this to work.
% 
% The optional 2nd output, if requested, is a struct with following fields:
%  st.FriendlyName (meaningful for Windows only)
%  st.Description  (OS dependent description)
%  st.MAC_address  (the same order as the 1st output)
%  st.IPv4_address (empty if not available)
%  st.IPv6_address (empty if not available)
% 
% Examples:
%  mac = MACAddress(); % return 1st MAC in string
% The format is like F0-4D-A2-DB-00-37 for Windows, f0:4d:a2:db:00:37 otherwise.
% 
%  macs = MACAddress(1); % return all MAC on the computer
% The output is cell even if only one MAC found.
% 
%  [macs, st] = MACAddress(1); % also return more info in st
% 
% To get numeric:
%  num = uint8(sscanf(MACAddress, '%2x%*c', 6))';

% 170510 Adapted this from RTBox code (Xiangrui.Li at gmail.com).
% 170525 Include mex for MS Windows.
% 171030 mex.c more robust. Include Octave 4 mex.
% 180525 use jsystem('ipconfig') for Windows; java method moved behind.
% 180626 implement 2nd optional output for both m and mex (almost rewritten).

if nargin<1 || isempty(allMac), allMac = false; end % default to first MAC

if ispc
    [tmp, str] = jsystem({'ipconfig.exe' '/all'});
	str = regexprep(str, '\r', '');
    mac_expr = 'Physical Address.*?:\s*((?:[0-9A-F]{2}-){5}[0-9A-F]{2})\s';
    nam_expr = '\nEthernet adapter\s+(.*?):?\n'; % nam/des/ip4/ip6 all in a block
    des_expr = 'Description.*?:\s*(.*?)\n';
    ip4_expr = 'IP(?:v4)? Address.*?:\s*((?:\d{1,3}\.){3}\d{1,3})';
    ip6_expr = 'IPv6 Address.*?:\s*((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    fmt = '%02X-%02X-%02X-%02X-%02X-%02X'; % adopt OS format preference
elseif ismac
    % [tmp, str] = jsystem({'networksetup' '-listallhardwareports'});
    [tmp, str] = jsystem({'ifconfig'});
    mac_expr = '\n\s+ether\s+((?:[0-9a-f]{2}:){5}[0-9a-f]{2})\s';
    des_expr = '\n(.*?):\s+';
    nam_expr = des_expr;
    ip4_expr = 'inet\s+((?:\d{1,3}\.){3}\d{1,3})';
    ip6_expr = 'inet6\s+((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    fmt = '%02x:%02x:%02x:%02x:%02x:%02x';
else % linux
    [err, str] = jsystem({'ip' 'address'}); % later Linux
    if ~err % almost always
        mac_expr = '\s+link/ether\s+((?:[0-9a-f]{2}:){5}[0-9a-f]{2})\s';
        des_expr = '\n\d+:\s+(.*?):\s+';
        ip4_expr = '\s+inet\s+((?:\d{1,3}\.){3}\d{1,3})';
        ip6_expr = '\s+inet6\s+((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    else % use ifconfig for old linux
        cmd = '/sbin/ifconfig';
        if ~exist(cmd, 'file'), cmd = 'ifconfig'; end
        [tmp, str] = jsystem({cmd});
        mac_expr = '\s+HWaddr\s+((?:[0-9a-f]{2}:){5}[0-9a-f]{2})\s';
        des_expr = '\n*(.*?)\s+';
        ip4_expr = 'inet addr:\s*((?:\d{1,3}\.){3}\d{1,3})';
        ip6_expr = 'inet6 addr:\s*((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    end
    nam_expr = des_expr;
    fmt = '%02x:%02x:%02x:%02x:%02x:%02x';
end

if allMac, [mac, ind] = regexp(str, mac_expr, 'tokens', 'start');
else,      [mac, ind] = regexp(str, mac_expr, 'tokens', 'start', 'once');
end
mac = [mac{:}];
% if iscell(mac) && numel(ind)>1 % make mac unique
%     [tmp, ia] = unique(mac);
%     ia = sort(ia); % [mac, ia] = unique(mac, 'stable'); ind = ind(ia);
%     mac = mac(ia);
%     ind = ind(ia);
% end

if nargout>1 && ~isempty(mac)
    st = struct('FriendlyName', [], 'Description', [], 'MAC_address', mac, ...
                'IPv4_address', [], 'IPv6_address', []);
    i0 = [1 regexp(str, '\n\S') numel(str)]; % split str into blocks
    for i = 1:numel(ind)
        j = find(i0<ind(i), 1, 'last');
        a = str(i0(j) : i0(j+1)); % the block with current mac
        c = regexp(a, nam_expr, 'tokens', 'once');
        if ~isempty(c), st(i).FriendlyName = c{1}; end
        c = regexp(a, des_expr, 'tokens', 'once');
        if ~isempty(c), st(i).Description  = c{1}; end
        c = regexp(a, ip4_expr, 'tokens', 'once');
        if ~isempty(c), st(i).IPv4_address = c{1}; end
        c = regexp(a, ip6_expr, 'tokens', 'once');
        if ~isempty(c), st(i).IPv6_address = c{1}; end
    end
end

% java method is OS-independent, more reliable than regexp, but often slower and
% miss eth1 seen at least 1 Ubuntu machine
if isempty(mac)
try %#ok
    if nargout>1, st = []; end
    if allMac, mac = {}; end
    ni = java.net.NetworkInterface.getNetworkInterfaces;
    while ni.hasMoreElements
        aa = ni.nextElement;
        a = aa.getHardwareAddress;
        if numel(a)~=6 || all(a==0), continue; end % not valid mac
        a = typecast(a, 'uint8'); % from int8
        m = sprintf(fmt, a);
        if nargout>1
            st(end+1).FriendlyName = char(aa.getName); %#ok
            st(end).Description = char(aa.getDisplayName);
            st(end).MAC_address = m;
            aa = aa.getInetAddresses;
            while aa.hasMoreElements
                c = char(aa.nextElement);
                a = regexp(c, '(\d{1,3}\.){3}\d{1,3}', 'match', 'once');
                if ~isempty(a), st(end).IPv4_address = a; end
                a = regexp(c, '(([0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})', 'match', 'once');
                if ~isempty(a)
                    st(end).IPv6_address = strrep(a, 'fe80:0:0:0:', 'fe80::'); 
                end
            end
        end
        if allMac, mac{end+1} = m; %#ok
        else, mac = m; break; % done after finding 1st
        end
    end
end
end

% If all attemps fail, give warning and return a random MAC
if isempty(mac)
    a = zeros(1,6);
    mac = sprintf(fmt, a);
    if nargout>1
        st = struct('FriendlyName', [], ...
            'Description', 'Failed to find network adapter', ...
            'MAC_address', mac, 'IPv4_address', [], 'IPv6_address', []);
    end
    if allMac, mac = {mac}; end
end

%% faster than system: based on https://github.com/avivrosenberg/matlab-jsystem
function [err, out] = jsystem(cmd)
% cmd is cell str, no quotation marks needed for file names with space.
try
    pb = java.lang.ProcessBuilder(cmd);
    pb.redirectErrorStream(true); % ErrorStream to InputStream
    process = pb.start();
    scanner = java.util.Scanner(process.getInputStream).useDelimiter('\A');
    if scanner.hasNext(), out = char(scanner.next()); else, out = ''; end
    err = process.exitValue; % err = process.waitFor() may hang
    if err, error('java.lang.ProcessBuilder error'); end
catch % fallback to system() if java fails like for Octave
    cmd = regexprep(cmd, '.+? .+', '"$0"'); % double quotes if with middle space
    [err, out] = system(sprintf('%s ', cmd{:}));
end

function [varargout] = judp(actionStr,varargin)
%
% judp.m--Uses Matlab's Java interface to handle User Datagram Protocol
% (UDP) communications with another application, either on the same
% computer or a remote one. 
%
% JUDP('SEND',PORT,HOST,MSSG) sends a message to the specifed port and
% host. HOST can be either a hostname (e.g., 'www.example.com') or a string
% representation of an IP address (e.g., '192.0.34.166'). Port is an
% integer port number between 1025 and 65535. The specified port must be
% open in the receiving machine's firewall.
% 
% MSSG = JUDP('RECEIVE',PORT,PACKETLENGTH) receives a message over the
% specified port. PACKETLENGTH should be set to the maximum expected
% message length to be received in the UDP packet; if too small, the
% message will be truncated.
%
% MSSG = JUDP('RECEIVE',PORT,PACKETLENGTH,TIMEOUT) attempts to receive a
% message but times out after TIMEOUT milliseconds. If TIMEOUT is not
% specified, as in the previous example, a default value is used.
%
% [MSSG,SOURCEHOST] = JUDP('RECEIVE',...) returns a string representation
% of the originating host's IP address, in addition to the received
% message. 
%
% Messages sent by judp.m are in int8 format. The int8.m function can be
% used to convert integers or characters to the correct format (use
% double.m or char.m to convert back after the datagram is received). 
% Non-integer values can be converted to int8 format using typecast.m (use
% typecast.m to convert back).
% 
% e.g.,   mssg = judp('receive',21566,10); char(mssg')
%         judp('send',21566,'208.77.188.166',int8('Howdy!'))         
%
% e.g.,   [mssg,sourceHost] = judp('receive',21566,10,5000)
%         judp('send',21566,'www.example.com',int8([1 2 3 4]))
%         
% e.g.,   mssg = judp('receive',21566,200); typecast(mssg,'double')
%         judp('send',21566,'localhost',typecast([1.1 2.2 3.3],'int8'))

% Developed in Matlab 7.8.0.347 (R2009a) on GLNX86.
% Kevin Bartlett (kpb@uvic.ca), 2009-06-18 16:11
%-------------------------------------------------------------------------

SEND = 1;
RECEIVE = 2;
DEFAULT_TIMEOUT = 1000; % [milliseconds]

% Handle input arguments.
if strcmpi(actionStr,'send')
    action = SEND;
    
    if nargin < 4
        error([mfilename '.m--SEND mode requires 4 input arguments.']);
    end % if
    
    port = varargin{1};
    host = varargin{2};
    mssg = varargin{3};
    
elseif strcmpi(actionStr,'receive')
    action = RECEIVE;
    
    if nargin < 3
        error([mfilename '.m--RECEIVE mode requires 3 input arguments.']);
    end % if
    
    port = varargin{1};
    packetLength = varargin{2};
    
    timeout = DEFAULT_TIMEOUT;
    
    if nargin > 3
        % Override default timeout if specified.
        timeout = varargin{3};
    end % if
    
else
    error([mfilename '.m--Unrecognised actionStr ''' actionStr ''.']);
end % if

% Test validity of input arguments.        
if ~isnumeric(port) || rem(port,1)~=0 || port < 1025 || port > 65535
    error([mfilename '.m--Port number must be an integer between 1025 and 65535.']);
end % if

if action == SEND
    if ~ischar(host)
        error([mfilename '.m--Host name/IP must be a string (e.g., ''www.examplecom'' or ''208.77.188.166''.).']);
    end % if
    
    if ~isa(mssg,'int8')
        error([mfilename '.m--Mssg must be int8 format.']);
    end % if
    
elseif action == RECEIVE    
    
    if ~isnumeric(packetLength) || rem(packetLength,1)~=0 || packetLength < 1
        error([mfilename '.m--packetLength must be a positive integer.']);
    end % if
    
    if ~isnumeric(timeout) || timeout <= 0
        error([mfilename '.m--timeout must be positive.']);
    end % if    
    
end % if

% Code borrowed from O'Reilly Learning Java, edition 2, chapter 12.
import java.io.*
import java.net.DatagramSocket
import java.net.DatagramPacket
import java.net.InetAddress

if action == SEND
    try
        addr = InetAddress.getByName(host);
        packet = DatagramPacket(mssg, length(mssg), addr, port);
        socket = DatagramSocket;
        socket.setReuseAddress(1);
        socket.send(packet);
        socket.close;
    catch 
        sendPacketError
        try
            socket.close;
        catch 
            closeError
            % do nothing.          
        end % try
        
        error('%s.m--Failed to send UDP packet.\nJava error message follows:\n%s',mfilename,sendPacketError.message);
        
    end % try
    
else
    try
        socket = DatagramSocket(port);
        socket.setSoTimeout(timeout);
        socket.setReuseAddress(1);
        packet = DatagramPacket(zeros(1,packetLength,'int8'),packetLength);        
        socket.receive(packet);
        socket.close;
        mssg = packet.getData;
        mssg = mssg(1:packet.getLength);     
        inetAddress = packet.getAddress;
        sourceHost = char(inetAddress.getHostAddress);
        varargout{1} = mssg;
        
        if nargout > 1
            varargout{2} = sourceHost;
        end % if
        
    catch 
        receiveError

        % Determine whether error occurred because of a timeout.
        if ~isempty(strfind(receiveError.message,'java.net.SocketTimeoutException'))
            errorStr = sprintf('%s.m--Failed to receive UDP packet; connection timed out.\n',mfilename);
        else
            errorStr = sprintf('%s.m--Failed to receive UDP packet.\nJava error message follows:\n%s',mfilename,receiveError.message);
        end % if        

        try
            socket.close;
        catch
            closeError
            % do nothing.
        end % try

        error(errorStr);
        
    end % try
    
end % if


