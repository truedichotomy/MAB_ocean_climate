function [err,L,U] = dg_ci(data)
    Z = 1.96; % for 0.95 CI
    stdev = std(data'); % std returns stdev of COLUMNS
    sem = stdev/sqrt(length(data(1,:)));
    L = - Z*sem';
    U = Z*sem';
    err = Z*sem';
