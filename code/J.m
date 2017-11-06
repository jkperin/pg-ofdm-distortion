function z = J(n,u)
% z = zeros(size(n));
% for k = 1:length(n)
%     if n(k) == 0
%         z(k) = 1 - u(k)^2/4;
%     elseif n(k) == 1 || n(k) == -1
%         z(k) = n(k)*u(k)/2;
%     else
%         z(k) = 0;
%     end
% end
z = besselj(n,u);
