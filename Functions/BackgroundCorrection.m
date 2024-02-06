function fit = BackgroundCorrection(flow,venc,order)
% Aaron Hess, university of oxford
% August 2017
% calculate the background phase correction map for 4D flow
% Determines which parts are static
% Fits a second order polynomial to the static flow data.
    time_dim = 4;
    flow_dim = 5;
    sz = size(flow);
    
    if(order>2)
        warning('BackgroundCorrection only tested for second order fit');
    end
    
    flow_mag = sqrt(sum(flow.^2,flow_dim));
    flow_std = std(flow_mag,[],time_dim);
    flow_none = and((flow_std<venc/30),(flow_std>0));  % SNR of 30 on average over all time fraemes
    flow_mean = mean(flow,time_dim);
    fit  =zeros([sz(1:3),1,sz(flow_dim)]);
    
    [X,Y,Z] = ind2sub(sz(1:3),find(flow_none));    
    ordStatic = ExpandOrd([X,Y,Z],order);
        
    [X, Y, Z] = ind2sub(size(flow_none),(1:prod(sz(1:3))).');
    ordAll = ExpandOrd([X,Y,Z],order);
    
    for iEnc = 1:sz(flow_dim)
        flow_mean_e = flow_mean(:,:,:,1,iEnc);        
        coef = ordStatic\flow_mean_e(flow_none==1);        
        fit(:,:,:,1,iEnc) = reshape(ordAll*coef, sz(1:3));
    end
end

function ord = ExpandOrd(ord1,order)
%     X2 = ord1(:,1).^2; Y2 = ord1(:,2).^2; Z2 = ord1(:,3).^2; XY = ord1(:,1).*ord1(:,2); XZ = ord1(:,1).*ord1(:,3); YZ = ord1(:,2).*ord1(:,3); C = ones(size(ord1,1),1);
%     ord = [C, ord1, X2, Y2, Z2, XY, XZ, YZ];

% Test code for nth order fit... not sure if it work
    ord = ord1;  % would be faster to pre-allocate ord.. just lazy
    for o = 2:order
        n = size(ord,2);
        for i = 1:n
            for j = i:n
                ord = [ord,ord(:,i).*ord(:,j)];
            end
        end
    end
    ord = [ord, ones(size(ord1,1),1)];
end