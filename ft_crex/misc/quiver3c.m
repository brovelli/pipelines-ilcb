function hq = quiver3c(pos, ori)
hq = quiver3(pos(:, 1), pos(:, 2), pos(:, 3), ori(:, 1), ori(:, 2), ori(:, 3));