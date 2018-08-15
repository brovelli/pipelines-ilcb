function wb_custcol(wb, col)
wbc = allchild(wb);
% Unlock java control for changing color
wbc(1).JavaPeer.setStringPainted(true);
wbc(1).JavaPeer.setForeground(java.awt.Color(col(1), col(2), col(3)));

set(wb.Children.Title, 'interpreter', 'none');