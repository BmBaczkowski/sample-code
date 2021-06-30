function sendMarker(dummymode, port, val)

if ~dummymode
pp_data(port, val);
WaitSecs(.002);
pp_data(port, 0);
end
end