function cleanup(dummymode_pp, pp)
  
  fclose('all');
  Eyelink('Shutdown');
  if ~dummymode_pp
  pp_close(pp);
  end
  sca;
  clear;
  ShowCursor;
  ListenChar(0);
  
end
