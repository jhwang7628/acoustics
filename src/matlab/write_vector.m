function write_vector (M, fname)
fid = fopen (fname, 'w');
n_row = size(M, 1);
fwrite (fid, n_row, 'int32');
fwrite (fid, M(:,1), 'double');
fclose (fid);
