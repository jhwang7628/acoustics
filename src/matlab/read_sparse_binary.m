function A = read_coord_sparse_matrix( fname )
    fid = fopen (fname, 'rb');
    n_row = fread (fid, 1, 'int32');
    
    if n_row == -1
        % Faster structs-of-arrays format
        n_row = fread (fid, 1, 'int32');
        n_col = fread (fid, 1, 'int32');
        n_nz = fread( fid, 1, 'int32' );
        rows = fread( fid, n_nz, 'int32' ) + 1;
        cols = fread( fid, n_nz, 'int32' ) + 1;
        vals = fread( fid, n_nz, 'double' );
        A = sparse( rows, cols, vals, n_row, n_col );        
    else
        % The slow array-of-structs format
        n_col = fread (fid, 1, 'int32');
        n_nz = fread( fid, 1, 'int32' );
        rows = zeros(n_nz,1);
        cols = zeros(n_nz,1);
        vals = zeros(n_nz,1);
        for i = 1:n_nz
            % NOTE: we need 1-based indices
            rows(i) = fread( fid, 1, 'int32' ) + 1;
            cols(i) = fread( fid, 1, 'int32' ) + 1; 
            vals(i) = fread( fid, 1, 'double' );
        end
        A = sparse( rows, cols, vals, n_row, n_col );
    end
    
    fclose( fid );
end
