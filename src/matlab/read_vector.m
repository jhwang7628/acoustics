function v = read_vector (fname)
fid = fopen (fname, 'r');
n_col = 1;
n_row = fread (fid, 1, 'int32');
v = fread (fid, [n_col n_row], 'double')';
fclose (fid);
