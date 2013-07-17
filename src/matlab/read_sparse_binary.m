function [ A, rows, cols, vals, n_row, n_col, n_nz ] = read_sparse_binary( ...
                                                          fname, adjust )
    fid = fopen (fname, 'rb');
    
    % Faster structs-of-arrays format
    n_row = fread (fid, 1, 'int32');
    n_col = fread (fid, 1, 'int32');
    n_nz = fread( fid, 1, 'int32' );
    rows = fread( fid, n_nz, 'int32' );
    cols = fread( fid, n_nz, 'int32' );
    vals = fread( fid, n_nz, 'double' );

    if ( nargin > 1 )
        rows = rows + 1;
        cols = cols + 1;
    end

    A = sparse( rows, cols, vals, n_row, n_col );        
    
    fclose( fid );
end
