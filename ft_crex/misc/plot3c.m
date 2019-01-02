function hp = plot3c(pos, varargin)
% Matlab plot3 function 
hpp = plot3(pos(:, 1), pos(:, 2), pos(:, 3), varargin{:});
if nargout
    hp = hpp;
end
