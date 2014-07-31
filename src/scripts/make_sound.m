function make_sound()
  p = read_vector('test.vector');

  % Normalize
  M = max(abs(p));
  p = p / (M * 1.01);

  % Using a high bitrate here since that is what's used for the newmark integrator
  % I'm pretty sure that Matlab provides a way of downsampling to a lower bitrate
  wavwrite(p, 192000, 'test.wav');
end
