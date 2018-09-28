function wb_custcol(wb, col)
wbc = allchild(wb);

% Unlock java control for changing color
wbc(1).JavaPeer.setStringPainted(true);
wbc(1).JavaPeer.setForeground(java.awt.Color(col(1), col(2), col(3)));
wbch = findobj(wb.Children, 'Type', 'Axes');
set(wbch(1).Title, 'interpreter', 'none');