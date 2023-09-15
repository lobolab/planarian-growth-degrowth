% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function [v] = hill(s, v_max, k_half, n) 
  v = (v_max * (s.^n)) ./ ((k_half)^n + (s).^n);
end
