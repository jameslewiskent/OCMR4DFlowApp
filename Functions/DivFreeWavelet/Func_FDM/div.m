function d = div(vx,vy,vz,res)
% d = div(vx,vy,vz,res)
% Calculates divergence

if (isempty(res))
    res = [1,1,1];
end
d = (vx-circshift(vx,[1,0,0]))/(res(1))+(vy-circshift(vy,[0,1,0]))/(res(2))+(vz-circshift(vz,[0,0,1]))/(res(3));
