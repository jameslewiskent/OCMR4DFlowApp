function out = Mozaic(img)
% mosaic 3D/ 4D imag

    sz = size(img);
    if(length(sz)<4)
        sz(4)=1;
    end
    n3D = sz(3);
    NbyN = ceil(sqrt(n3D));
    out= zeros(sz(1)*NbyN,sz(2)*NbyN,sz(4));

    start_row = 1:sz(1):sz(1)*NbyN;
    end_row = start_row + sz(1)-1;
    start_col = 1:sz(2):sz(2)*NbyN;
    end_col = start_col + sz(2)-1;

    for i = 1:sz(3)
        icol = mod(i-1,NbyN)+1;
        irow = ceil(i/NbyN);
        out(start_row(irow):end_row(irow),start_col(icol):end_col(icol),:) = squeeze(abs(img(:,:,i,:)));
    end
end