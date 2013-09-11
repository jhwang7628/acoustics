function make_material_parameters(youngsModulus, poissonRatio, filename)
  materialParameters = [ youngsModulus; poissonRatio ];

  write_vector( materialParameters, filename );
end
