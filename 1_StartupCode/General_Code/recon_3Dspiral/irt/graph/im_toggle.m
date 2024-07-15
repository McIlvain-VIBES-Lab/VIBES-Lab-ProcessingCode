  function im_toggle(i1, i2, varargin)
%|function im_toggle(i1, i2, [..., options for im()])
%| toggle between two or more images via keypress

if nargin == 1 && streq(i1, 'test'), im_toggle_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% find leading additional arguments corresponding to images
iall = {i1, i2};
names = {inputname(1), inputname(2)};
ii = 3;
while length(varargin) && isequal(size(i2), size(varargin{1}))
	iall = {iall{:}, varargin{1}};
	varargin = {varargin{2:end}};
	names{ii} = inputname(ii);
	ii = ii + 1;
end

if ~im, return, end

% toggle between two or more images

ft.args = {{}, {'horiz', 'right'}};
ft.pos = [0.01, 0.99];

while (1)
	im clf

	for ii=1:length(iall)
		im clf
		im(iall{ii}, varargin{:})
		fig_text(ft.pos(2-mod(ii,2)), 0.01, ...
			sprintf(['toggle i%d: ' inputname(ii)], ii), ...
			ft.args{2-mod(ii,2)})

%		pause
		in = input('hit enter for next image, or "q" to quit ', 's');
		if streq(in, 'q')
			set(gca, 'nextplot', 'replace')
			return
		end
	end
end


function im_toggle_test
nx = 20;
i1 = eye(20);
i2 = flipud(i1);
i3 = 1 - i1;
im_toggle(i1, i2, i3, [0 2])
