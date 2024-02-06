function [ ims_comb , rx_maps ] = rx_combination( ims )
% RX_COMBINATION Optimally combines multi-rx data. Input ims is an Nd
% array where the first three dims are three spatial ordinates, and the
% last dimension is the rx coil dimension. All dimension between the 3rd
% and last will be used to estimate Rx sensitivities. The Brunner SVD approach
% is used here (ISMRM2010 #272).

sz = size(ims);

rx_maps = zeros([sz(1:3) sz(end)]);
ims_comb = zeros(sz(1:(end-1)));
for x = 1:sz(1)
    for y = 1:sz(2)
        for z = 1:sz(3)
            tmp = reshape(ims(x,y,z,:),[prod(sz(4:(end-1))) sz(end)]);
            [~,~,v] = svd(tmp);
            rx_maps(x,y,z,:) = permute(v(:,1),[2 3 4 1]);
            ims_comb(x,y,z,:) = tmp*v(:,1);
        end
    end
end
ims_comb = reshape(ims_comb,sz(1:(end-1)));
end

