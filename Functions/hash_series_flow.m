function hash = hash_series_flow(ser_flow)
% Hash series flow IDs for filename
hasher = System.Security.Cryptography.HashAlgorithm.Create('MD5');
hash_byte = hasher.ComputeHash( uint8(ser_flow) );  % System.Byte class
hash_uint8 = uint8( hash_byte );               % Array of uint8
hash_hex = dec2hex(hash_uint8);                % Array of 2-char hex codes

% Generate the hex codes as 1 long series of characters
hashStr = string([]);
nBytes = length(hash_hex);
for k=1:nBytes
    hashStr(end+1:end+2) = hash_hex(k,:);
end
hash = char(strjoin(hashStr,''));

end

