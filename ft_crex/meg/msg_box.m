function msg_box(smsg, tit)
if nargin < 2
    tit = 'Message';
end

sz = 12;

cmsg = char(smsg);
nc = length(cmsg(1,:));
wbox = round(nc*sz/1.7);
nr = length(smsg);

hb = msgbox(repmat({' '}, nr+3, 1), tit);
hb.Color = [1 1 1];
pos = hb.Position;

hb.Position = [pos(1)-round(wbox/4) pos(2)-100 wbox pos(4)];

hok = findobj(hb, 'Tag', 'OKButton');
hok.BackgroundColor = [1 0.98 0.96];
op = hok.Position;
hok.Position = [round(wbox/2 - op(3)/2) op(2:4)];
htxt = findobj(hb, 'Tag', 'MessageBox');
htxt.FontSize = 12;
htxt.VerticalAlignment = 'bottom';
htxt.String = [smsg ; ' '];