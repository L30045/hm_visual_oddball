% save multiple copies of figure in different types
% inputs:
%     path
%     filetype or cell of filetypes
% ex: print_multi('~/figures', {'png', 'pdf'})

function print_multi(varargin)
if nargin < 2 || ischar(varargin{2})
    print(varargin{:});
else
    % check inputs
    types = varargin{2};
    assert(iscellstr(types), ...
        'Third input must be either string or cell string')
%     if any(strcmp(types, 'eps')) && ~any(types, 'epsc')
%         warning('saveas_multi: are you sure you want ''eps'' and not ''epsc''? ''epsc'' produces color images while ''eps'' produces black and white images');
%     end
    % set paper position
%     set(gcf, 'paperunits', 'points')
%     set(gcf, 'paperposition', get(gcf, 'position'))
    set(gcf, 'PaperPositionMode', 'auto')
    % print requested types
    for it = 1:length(types)
        if strcmp(types{it}, 'pdf')
            print(varargin{1}, '-dpdf', '-r0')
        elseif strncmp(types{it}, 'eps', 3)
            print(varargin{1}, '-depsc2', '-r0')
        else
            print(varargin{1}, ['-d' types{it}])
        end
    end
end