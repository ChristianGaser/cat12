% CPRINTF displays styled formatted text in the Command Window
%
% Syntax:
%    count = cprintf(style,format,...)
%
% Description:
%    CPRINTF processes the specified text using the exact same FORMAT
%    arguments accepted by the built-in SPRINTF and FPRINTF functions.
%
%    CPRINTF then displays the text in the Command Window using the
%    specified STYLE argument. The accepted styles are those used for
%    Matlab's syntax highlighting (see: File / Preferences / Colors / 
%    M-file Syntax Highlighting Colors), and also user-defined colors.
%
%    The possible pre-defined STYLE names are (case-insensitive):
%
%      'Text'                - default (based on Preferences): black
%      'Keywords'            - default (based on Preferences): blue
%      'Comments'            - default (based on Preferences): green
%      'Strings'             - default (based on Preferences): purple
%      'UnterminatedStrings' - default (based on Preferences): dark red
%      'SystemCommands'      - default (based on Preferences): orange
%      'Errors'              - default (based on Preferences): light red
%      'Hyperlinks'          - default (based on Preferences): blue (underlined)
%
%      'Black','Cyan','Magenta','Blue','Green','Red','Yellow','White',
%      'Gray','Orange','Pink','Silver','Maroon','Purple','Olive','Navy',
%      'Teal','Violet','Rose','Brown','Indigo','Gold','Wheat'
%       ...and also variants of these colors using a 'light' or 'dark' prefix.
%       For example: 'light-gray', 'darkGreen'
%
%    STYLE is case-insensitive and accepts unique (non-ambiguous) partial strings
%       For example: 'cy' [instead of 'cyan'], 'system', 'err'.
%
%    STYLE custom colors can be specified in 3 variants:
%       [0.1,0.7,0.3] - standard Matlab RGB color format (range: 0.0-1.0)
%       [26, 178, 76] - numeric RGB values (range: 0-255)
%       '#1ab34d'     - Hexadecimal format (range: '00'-'FF', case insensitive)
%                       3-digit HTML format is also accepted: 'a5f'='aa55ff'
%
%    STYLE beginning with '-' or '_' will be underlined.
%       For example: '-Blue', '_Comments', '-#0FF', -[0,1,1]
%
%    STYLE beginning with '*' will be bold (on R2011b+ only).
%       For example: '*Blue', '*Comments', '*#F00', '*[1,0,0]'
%       Bold & underline can be combined, for example: '*_red' or '_*blue'
%
%    CPRINTF by itself, without any input parameters, displays the usage demo
%
% Usage examples:
%    cprintf;   % displays the usage demo
%    cprintf('text',       'regular black text');
%    cprintf('hyper',      'followed %s','by');
%    cprintf('key',        '%d colored', 8);
%    cprintf('-comment',   'and underlined');
%    cprintf('err',        'elements:\n');
%    cprintf('cyan',       'cyan ');
%    cprintf('_green',     'underlined green ');
%    cprintf('lightgreen', 'light green ');
%    cprintf('dark-green', 'dark green ');
%    cprintf('#800080',    'any custom RGB ');
%    cprintf([1,0.5,0],    'multi-\nline orange ');
%    cprintf(-[1,0,1],     'underlined magenta\n');
%    cprintf('*blue',      'and *bold* (R2011b+) ');
%    cprintf('*_red',      'and *bold underline* (R2011b+)\n');
%
%    cprintf('string');  % same as fprintf('string') and cprintf('text','string')
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab
%    functionality. It works on Matlab 7+, but use at your own risk!
%
%    A technical description of the implementation can be found at:
%    <a href="http://undocumentedmatlab.com/articles/cprintf">http://UndocumentedMatlab.com/articles/cprintf</a>
%
% Limitations:
%    1. In R2011a and earlier, a single space char is inserted at the
%       beginning of each CPRINTF text segment (this is ok in R2011b+).
%
%    2. In R2011a and earlier, consecutive differently-colored multi-line
%       CPRINTFs sometimes display incorrectly on the bottom line.
%       As far as I could tell this is due to a Matlab bug. Examples:
%         >> cprintf('-str','under\nline'); cprintf('err','red\n'); % hidden 'red', non-hidden '_'
%         >> cprintf('str','regu\nlar');    cprintf('err','red\n'); % underline red (not purple) 'lar'
%
%    3. Sometimes, non newline ('\n')-terminated segments display unstyled
%       (black) when the command prompt chevron ('>>') regains focus on the
%       continuation of that line (I can't pinpoint when this happens). 
%       To fix this, simply newline-terminate all command-prompt messages.
%
%    4. In R2011b and later, the above errors appear to be fixed. However,
%       the last character of an underlined segment is not underlined for
%       some unknown reason (add an extra space character to make it look better)
%
%    5. In old Matlab versions (e.g., Matlab 7.1 R14), multi-line styles
%       only affect the first line. Single-line styles work as expected.
%       R14 also appends a single space after underlined segments.
%
%    6. Bold style is only supported on R2011b+ (can be combined with underline).
%
%    7. CPRINTF is not supported in Matlab mobile, Live Editor, deployed, diary,
%       and terminal (no desktop) modes; The new web-based (JavaScript) desktop
%       is only supported in R2025a or newer. This limitation depends on internal
%       Matlab limitations that may possibly be lifted in future Matlab releases.
%
%    8. On R2025a or newer, HTML hyperlinks cannot have a non-default style.
%
% See also:
%    sprintf, fprintf

% Change log:
%    2009-05-13: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>
%    2009-05-28: corrected nargout behavior suggested by Andreas Gäb
%    2009-09-28: Fixed edge-case problem reported by Swagat K
%    2010-06-27: Fix for R2010a/b; fixed edge case reported by Sharron; CPRINTF with no args runs the demo
%    2011-03-04: Performance improvement
%    2011-08-29: Fix by Danilo (FEX comment) for non-default text colors
%    2011-11-27: Fixes for R2011b
%    2012-08-06: Fixes for R2012b; added bold style; accept RGB string (non-numeric) style
%    2012-08-09: Graceful degradation support for deployed (compiled) and non-desktop applications; minor bug fixes
%    2015-03-20: Fix: if command window isn't defined yet (startup) use standard fprintf as suggested by John Marozas
%    2015-06-24: Fixed a few discoloration issues (some other issues still remain)
%    2020-01-20: Fix by T. Hosman for embedded hyperlinks
%    2021-04-07: Enabled specifying color as #RGB (hexa codes), [.1,.7,.3], [26,178,76]
%    2022-01-04: Fixed cases of invalid colors (especially bad on R2021b onward)
%    2022-03-26: Fixed cases of using string (not char) inputs
%    2025-02-12: Support R2025a (web-based desktop)
%    2025-03-05: Output to STDERR if style is 'error','red' or [1,0,0] (non-Desktop modes only)
%    2025-09-07: Fixed auto-hyperlink in R2025a+; empty style now means 'text'
%    2025-09-08: Fixed "stuck" color styling in case of invalid escape sequence reported by Alex N
%    2026-02-02: Workaround for hyperlink style on Javascript & online desktops reported by Meike
%    2026-04-20: Added predefined color names & Dark/Light variants, updated demo
%    2026-04-22: Added several predefined named colors; reverted unnecessary R2026a fix
%    2026-06-26: Added ability for concurrent bold & underline (thanks @SongYuxuan)

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.22 $  $Date: 2026/06/26 14:24:00 $
