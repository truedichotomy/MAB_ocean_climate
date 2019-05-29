% DG 2017-08-30, 2017-11-01, 2017-12-15, 2018-01-04
% this script perform regression analysis for time series data from the MAB taking potential autocorrelation (effective DOF) into account.
% need to run dg_grid_regions_hydroavg.m before this script to calculate the average temperature and salinity for each region of each season.

% setup directories needed for the script
dg_setup_MABclimate_dir

loadmonthlyflag = 1
loadhydroflag = 1
calcflag = 1
seasonbiasflag = 1 % 0 (without seasonalbias) or 1 (with seasonal bias correction)
dataflag = 'depth-averaged' % 'bottom' or 'depth-averaged' or 'bottomSODA' or 'depth-averagedSODA', only the first and last letters matter
saveflag = 1
sodaflag = 0

searchradius = 40000
p = 4

% variable (4): 1-temp, 2-tempshw, 3-salt, 4-saltshw
% region (9): 1-SNE, 2-NYB1, 3-NYB2, 4-SS1, 5-SS2, 6-MAB, 7-GB, 8-ENE, 9-GOM
% season (3): 1-Winter-Spring, 2-Spring-Summer, 3-Fall-Winter
% time period (3): 1-1977-1999, 2-1999-2017, 3-1977-2016

%si = 1 % season (3): 1-Winter-Spring, 2-Spring-Summer, 3-Fall-Winter
%ri = 1 % region (9): 1-SNE, 2-NYB1, 3-NYB2, 4-SS1, 5-SS2, 6-MAB, 7-GB, 8-ENE, 9-GOM
%vi = 1 % variable (4): 1-temp, 2-tempshw, 3-salt, 4-saltshw
%ti = 3 % time period (3): 1-1977-1999, 2-1999-2017, 3-1977-2016

ns = 3;
nr = 9;
nv = 4;
nt = 3;

if sodaflag == 0
	tind{1} = [1:23];
	tind{2} = [23:40];
	tind{3} = [1:40];
elseif sodaflag == 1
	tind{1} = [4:23];
	tind{2} = [26:37];
	tind{3} = [4:39];
end %if

switch seasonbiasflag
case 0
	var{1} = 'temp';
	var{2} = 'tempshw';
	var{3} = 'salt';
	var{4} = 'saltshw';
case 1
	var{1} = 'temp_nobias';
	var{2} = 'tempshw_nobias';
	var{3} = 'salt_nobias';
	var{4} = 'saltshw_nobias';
end %switch

m = repmat(NaN,[ns,nr,nv,nt]);
b = m;
mloeff = m;
mhieff = m;
mlo = m;
mhi = m;
tsN = m;
tsNstar = m;
tsNeff = m;
hypothesis = m;
lseason = m;
lregion = m;
lvariable = m;
ltimeperiod = m;

if loadmonthlyflag == 1
	if ~exist('casts')
		load([workdir 'hydroMABcasts.mat']);
	end %if

	if ~exist('region')
		region = dg_grid_regions_define;
	end %if

	if sodaflag == 0
		hydroMABavgmonthlypath = [workdirlocal 'hydroMABavgmonthly_r' num2str(searchradius) '_p' num2str(p) '.mat'];
	elseif sodaflag == 1
		hydroMABavgmonthlypath = [workdirlocal 'hydroMABavgmonthlySODA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
	end %if
	load(hydroMABavgmonthlypath);

	yyyylist = [1977:2016];
	month = [1:12];
	regions = [1:9];

	mtemp = repmat(NaN,[length(month),length(regions)]);
	msalt = mtemp;
	mtempshw = mtemp;
	msaltshw = mtemp;
	mvoltemp = mtemp;
	mvolsalt = mtemp;
	mvolshw = mtemp;
	mvolume = mtemp;

	% load regional boundary data and find the time and location of each cast
	reg = dg_grid_regions_define_basic;
	%region = dg_grid_regions_define;
	yyyy = [casts.yr]';
	dyd = [casts.dyd]';
	yd = [casts.yd]';
	mm = str2num(datestr(datenum(yyyy,1,dyd),5));
	londata = [casts.lon]';
	latdata = [casts.lat]';
	if sodaflag == 1
		lonind = find(londata > 180);
		londata(lonind) = londata(lonind) - 360;
	end %if

	depth = -double([casts.zdeep])';

	for ri = 1:9
		for mi = 1:12
			mtemp(mi,ri) = hydromom.month(mi).region(ri).mtemp;
			msalt(mi,ri) = hydromom.month(mi).region(ri).msalt;
			mtempshw(mi,ri) = hydromom.month(mi).region(ri).mtempshw;
			msaltshw(mi,ri) = hydromom.month(mi).region(ri).msaltshw;
			mvoltemp(mi,ri) = hydromom.month(mi).region(ri).voltemp;
			mvolsalt(mi,ri) = hydromom.month(mi).region(ri).volsalt;
			mvolume(mi,ri) = hydromom.month(mi).region(ri).volume;
		end %for mi
	end %for ri

	ccc.month = month;
	ccc.tempmonthly = mtemp;
	ccc.saltmonthly = msalt;
	ccc.tempshwmonthly = mtempshw;
	ccc.saltshwmonthly = msaltshw;
	ccc.voltempmonthly = mvoltemp;
	ccc.volsaltmonthly = mvolsalt;
	ccc.volumemonthly = mvolume;
    ccc.ncastsmonthly = repmat(NaN,[length(ccc.month), length(reg)]);
	ccc.ncasts = repmat(NaN,[length(yyyylist), length(ccc.month), length(reg)]);

	for ri = 1:length(reg)
	    inregion = inpolygon(londata,latdata,reg(ri).lon,reg(ri).lat);
	    regionalind = find(inregion == logical(1));
	    reg(ri).clon = londata(regionalind);
	    reg(ri).clat = latdata(regionalind);
	    reg(ri).cdepth = depth(regionalind);
		reg(ri).cyyyy = yyyy(regionalind);
	    reg(ri).cmm = mm(regionalind);
	    reg(ri).ind = [];

	    for mi = 1:12
	        %reg(ri).ind{mi} = find(cmm == mi & cdepth <= 1500);
	        ccc.ncastsmonthly(mi,ri) = length(find(reg(ri).cmm == mi & reg(ri).cdepth <= 1500));
	    end %for

		for yi = 1:length(yyyylist)
			for mi = 1:12
				%reg(ri).ind{mi} = find(cmm == mi & cdepth <= 1500);
				ccc.ncasts(yi,mi,ri) = length(find(reg(ri).cmm == mi & reg(ri).cyyyy == yyyylist(yi) & reg(ri).cdepth <= 1500));
			end %for
		end %for
	end %for

	if saveflag == 1
		if sodaflag == 0
			save([workdirlocal 'MABmonthly2D_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydromom','ccc','reg','month','yyyylist','-v7.3')
		elseif sodaflag == 1
			save([workdirlocal 'MABmonthly2D_SODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydromom','ccc','reg','month','yyyylist','-v7.3')
		end %if
	end %if

end %if loadmonthlyflag

if loadhydroflag == 1
	if lower(dataflag(1)) == 'd' & dataflag(end) ~= 'A'
		hydroMABavgpath = [workdirlocal 'hydroMABavg2Da_r' num2str(searchradius) '_p' num2str(p) '.mat'];
		load(hydroMABavgpath);
	elseif lower(dataflag(1)) == 'b' & dataflag(end) ~= 'A'
		hydroMABavgpath = [workdirlocal 'hydroMABavg2Db_r' num2str(searchradius) '_p' num2str(p) '.mat'];
		load(hydroMABavgpath);
		hydrom = hydrob;
	elseif lower(dataflag(1)) == 'd' & dataflag(end) == 'A'
		hydroMABavgpath = [workdirlocal 'hydroMABavg2DaSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
		load(hydroMABavgpath);
	elseif lower(dataflag(1)) == 'b' & dataflag(end) == 'A'
		hydroMABavgpath = [workdirlocal 'hydroMABavg2DbSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'];
		load(hydroMABavgpath);
		hydrom = hydrob;
	end %if

	if seasonbiasflag == 1
		%nobiasMABpath = [workdir 'nobias_ts_2D.mat'];
		%nbccc = load(nobiasMABpath);
	end %if

	if ~exist('region')
		region = dg_grid_regions_define;
	end %if

	yyyy = [1977:2016]';
	ccc.yyyy = yyyy;

	%season(1).label = 'Year Day 1-125';
	%season(2).label = 'Year Day 126-250';
	%season(3).label = 'Year Day 251-365';
	season(1).label = 'Year Day 1-120';
	season(2).label = 'Year Day 121-243';
	season(3).label = 'Year Day 244-365';
	season(1).name = 'Winter-Spring';
	season(2).name = 'Spring-Summer';
	season(3).name = 'Fall-Winter';
	season(1).month = [1:4];
	season(2).month = [5:8];
	season(3).month = [9:12];

	ccc.season = season;
	ccc.temp = repmat(NaN,[length(yyyy), 3, 9]);
	ccc.salt = repmat(NaN,[length(yyyy), 3, 9]);
	ccc.tempshw = ccc.temp;
	ccc.saltshw = ccc.salt;
	ccc.voltemp = ccc.temp;
	ccc.volsalt = ccc.salt;
	ccc.volshw = ccc.salt;
	ccc.tempfilt = ccc.temp;
	ccc.saltfilt = ccc.salt;
	ccc.tempshwfilt = ccc.temp;
	ccc.saltshwfilt = ccc.salt;
	ccc.tempavg = repmat(NaN,[3, 9]);
	ccc.tempstd = ccc.tempavg;
	ccc.saltavg = ccc.tempavg;
	ccc.saltstd = ccc.tempavg;
	ccc.tempshwavg = ccc.tempavg;
	ccc.tempshwstd = ccc.tempavg;
	ccc.saltshwavg = ccc.tempavg;
	ccc.saltshwstd = ccc.tempavg;

	ccc.tempfilt_nobias = repmat(NaN,[length(yyyy), 3, 9]);;
	ccc.saltfilt_nobias = repmat(NaN,[length(yyyy), 3, 9]);;
	ccc.tempshwfilt_nobias = repmat(NaN,[length(yyyy), 3, 9]);;
	ccc.saltshwfilt_nobias = repmat(NaN,[length(yyyy), 3, 9]);;
	ccc.tempavg_nobias = repmat(NaN,[3, 9]);
	ccc.tempstd_nobias = repmat(NaN,[3, 9]);
	ccc.saltavg_nobias = repmat(NaN,[3, 9]);
	ccc.saltstd_nobias = repmat(NaN,[3, 9]);
	ccc.tempshwavg_nobias = repmat(NaN,[3, 9]);
	ccc.tempshwstd_nobias = repmat(NaN,[3, 9]);
	ccc.saltshwavg_nobias = repmat(NaN,[3, 9]);
	ccc.saltshwstd_nobias = repmat(NaN,[3, 9]);


	for yi = 1:length(yyyy)
	    for si = 1:3
	        [yi si]
	        ccc.temp(yi,si,:) = [hydrom(yi).season(si).region(:).temp];
	        ccc.salt(yi,si,:) = [hydrom(yi).season(si).region(:).salt];
	        ccc.tempshw(yi,si,:) = [hydrom(yi).season(si).region(:).tempshw];
	        ccc.saltshw(yi,si,:) = [hydrom(yi).season(si).region(:).saltshw];
	        ccc.voltemp(yi,si,:) = [hydrom(yi).season(si).region(:).voltemp];
	        ccc.volsalt(yi,si,:) = [hydrom(yi).season(si).region(:).volsalt];
	        ccc.volshw(yi,si,:) = [hydrom(yi).season(si).region(:).volshw];
	    end %for si
	end %for yi

	for si = 1:ns
	    for ri = 1:nr
	        ccc.tempfilt(:,si,ri) = dg_extrema_remove(ccc.temp(:,si,ri),1);
	        ccc.saltfilt(:,si,ri) = dg_extrema_remove(ccc.salt(:,si,ri),1);
	        ccc.tempshwfilt(:,si,ri) = dg_extrema_remove(ccc.tempshw(:,si,ri),1);
	        ccc.saltshwfilt(:,si,ri) = dg_extrema_remove(ccc.saltshw(:,si,ri),1);

	        ccc.tempavg(si,ri) = nanmean(ccc.tempfilt(:,si,ri));
	        ccc.tempstd(si,ri) = nanstd(ccc.tempfilt(:,si,ri));
	        ccc.saltavg(si,ri) = nanmean(ccc.saltfilt(:,si,ri));
	        ccc.saltstd(si,ri) = nanstd(ccc.saltfilt(:,si,ri));
	        ccc.tempshwavg(si,ri) = nanmean(ccc.tempshwfilt(:,si,ri));
	        ccc.tempshwstd(si,ri) = nanstd(ccc.tempshwfilt(:,si,ri));
	        ccc.saltshwavg(si,ri) = nanmean(ccc.saltshwfilt(:,si,ri));
	        ccc.saltshwstd(si,ri) = nanstd(ccc.saltshwfilt(:,si,ri));
	    end %for
	end %for

	% seasonal bias correction calculations
	ccc.deltaT = repmat(NaN,[length(yyyy),ns,nr]);
	ccc.deltaS = repmat(NaN,[length(yyyy),ns,nr]);
	ccc.deltaTshw = repmat(NaN,[length(yyyy),ns,nr]);
	ccc.deltaSshw = repmat(NaN,[length(yyyy),ns,nr]);
	for yi = 1:length(ccc.yyyy)
		for si = 1:ns
			for ri = 1:nr
				ccc.deltaT(yi,si,ri) = nanmean(ccc.tempmonthly(season(si).month,ri)) - nansum(ccc.tempmonthly(season(si).month,ri) .* ccc.ncasts(yi,season(si).month,ri)') / nansum(ccc.ncasts(yi,season(si).month,ri));

				ccc.deltaS(yi,si,ri) = nanmean(ccc.saltmonthly(season(si).month,ri)) - nansum(ccc.saltmonthly(season(si).month,ri) .* ccc.ncasts(yi,season(si).month,ri)') / nansum(ccc.ncasts(yi,season(si).month,ri));

				ccc.deltaTshw(yi,si,ri) = nanmean(ccc.tempshwmonthly(season(si).month,ri)) - nansum(ccc.tempshwmonthly(season(si).month,ri) .* ccc.ncasts(yi,season(si).month,ri)') / nansum(ccc.ncasts(yi,season(si).month,ri));

				ccc.deltaSshw(yi,si,ri) = nanmean(ccc.saltshwmonthly(season(si).month,ri)) - nansum(ccc.saltshwmonthly(season(si).month,ri) .* ccc.ncasts(yi,season(si).month,ri)') / nansum(ccc.ncasts(yi,season(si).month,ri));
			end %for
		end %for
	end %for


	if seasonbiasflag == 1
		%cccLW.temp_nobias = nbccc.ccc.temp_nobias;
		%cccLW.tempshw_nobias = nbccc.ccc.tempshw_nobias;
		%cccLW.salt_nobias = nbccc.ccc.salt_nobias;
		%cccLW.saltshw_nobias = nbccc.ccc.saltshw_nobias;

		ccc.temp_nobias = ccc.temp + ccc.deltaT;
		ccc.tempshw_nobias = ccc.tempshw + ccc.deltaTshw;
		ccc.salt_nobias = ccc.salt + ccc.deltaS;
		ccc.saltshw_nobias = ccc.saltshw + ccc.deltaSshw;

		for si = 1:ns
			for ri = 1:nr
				ccc.tempfilt_nobias(:,si,ri) = dg_extrema_remove(ccc.temp(:,si,ri),1);
				ccc.saltfilt_nobias(:,si,ri) = dg_extrema_remove(ccc.salt(:,si,ri),1);
				ccc.tempshwfilt_nobias(:,si,ri) = dg_extrema_remove(ccc.tempshw(:,si,ri),1);
				ccc.saltshwfilt_nobias(:,si,ri) = dg_extrema_remove(ccc.saltshw(:,si,ri),1);

				ccc.tempavg_nobias(si,ri) = nanmean(ccc.tempfilt_nobias(:,si,ri));
				ccc.tempstd_nobias(si,ri) = nanstd(ccc.tempfilt_nobias(:,si,ri));
				ccc.saltavg_nobias(si,ri) = nanmean(ccc.saltfilt_nobias(:,si,ri));
				ccc.saltstd_nobias(si,ri) = nanstd(ccc.saltfilt_nobias(:,si,ri));
				ccc.tempshwavg_nobias(si,ri) = nanmean(ccc.tempshwfilt_nobias(:,si,ri));
				ccc.tempshwstd_nobias(si,ri) = nanstd(ccc.tempshwfilt_nobias(:,si,ri));
				ccc.saltshwavg_nobias(si,ri) = nanmean(ccc.saltshwfilt_nobias(:,si,ri));
				ccc.saltshwstd_nobias(si,ri) = nanstd(ccc.saltshwfilt_nobias(:,si,ri));
			end %for
		end %for

	end %if

	if lower(dataflag(1)) == 'd' & dataflag(end) ~= 'A'
		ccc.dataflag = 'depth-averaged'
		switch seasonbiasflag
		case 0
		    save([workdirlocal 'MABclimate2Da_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
		case 1
			save([workdirlocal 'MABclimate2Da_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
		end %switch
	elseif lower(dataflag(1)) == 'b' & dataflag(end) ~= 'A'
		ccc.dataflag = 'bottom'
		switch seasonbiasflag
		case 0
		    save([workdirlocal 'MABclimate2Db_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
		case 1
			save([workdirlocal 'MABclimate2Db_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
		end %switch
	elseif lower(dataflag(1)) == 'd' & dataflag(end) == 'A'
		ccc.dataflag = 'depth-averaged-SODA'
		switch seasonbiasflag
		case 0
		    save([workdirlocal 'MABclimate2DaSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
		case 1
			save([workdirlocal 'MABclimate2DaSODA_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
		end %switch
	elseif lower(dataflag(1)) == 'b' & dataflag(end) == 'A'
		ccc.dataflag = 'bottom-SODA'
		switch seasonbiasflag
		case 0
		    save([workdirlocal 'MABclimate2DbSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
		case 1
			save([workdirlocal 'MABclimate2DbSODA_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
		end %switch
	end %if
	display('Done loading hydro data.')
end %if loadhydroflag

if calcflag == 1
	for si = 1:ns
	    for ri = 1:nr
	        [season(si).name '  ' region(ri).label]
	        for vi = 1:nv
	            for ti = 1:nt
					%si = 3; ri = 9; vi = 3; ti = 2;
	                x = ccc.yyyy(tind{ti});
	                y = eval(['ccc.' var{vi} '(tind{ti},si,ri)']);
					%badind = find((detrend(y) > 3.2*nanstd(y)) | (detrend(y) < -3.2*nanstd(y)));
					%if ~isempty(badind)
					%	y(badind) = NaN;
					%end %if
	                rega = dg_DFanalysis_core(x, y);
					b(si,ri,vi,ti) = rega.b;
	                m(si,ri,vi,ti) = rega.m;
					mlo(si,ri,vi,ti) = rega.mlo;
	                mhi(si,ri,vi,ti) = rega.mhi;
	                mloeff(si,ri,vi,ti) = rega.mloeff;
	                mhieff(si,ri,vi,ti) = rega.mhieff;
					tsN(si,ri,vi,ti) = rega.N;
					tsNstar(si,ri,vi,ti) = rega.Nstar;
					tsNeff(si,ri,vi,ti) = rega.Neff;
	                hypothesis(si,ri,vi,ti) = rega.hypothesis;
					lseason(si,ri,vi,ti) = si;
					lregion(si,ri,vi,ti) = ri;
					lvariable(si,ri,vi,ti) = vi;
					ltimeperiod(si,ri,vi,ti) = ti;
	            end %for ti
	        end %for vi
	    end %for ri
	end %for si

	ccc.mlabel = 'season,region,variable,timeperiod';
	ccc.season = lseason;
	ccc.region = lregion;
	ccc.variable = lvariable;
	ccc.timeperiod = ltimeperiod;
	ccc.b = b;
	ccc.m = m;
	ccc.mlo = mlo;
	ccc.mhi = mhi;
	ccc.mloeff = mloeff;
	ccc.mhieff = mhieff;
	ccc.N = tsN;
	ccc.Nstar = tsNstar;
	ccc.Neff = round(tsNeff);
	ccc.hypothesis = hypothesis;

	ccc.season1d = reshape(lseason,[],1);
	ccc.region1d = reshape(lregion,[],1);
	ccc.variable1d = reshape(lvariable,[],1);
	ccc.timeperiod1d = reshape(ltimeperiod,[],1);

	ccc.b1d = reshape(b,[],1);
	ccc.m1d = reshape(m,[],1);
	ccc.mlo1d = reshape(mlo,[],1);
	ccc.mhi1d = reshape(mhi,[],1);
	ccc.mloeff1d = reshape(mloeff,[],1);
	ccc.mhieff1d = reshape(mhieff,[],1);
	ccc.N1d = reshape(tsN,[],1);
	ccc.Nstar1d = reshape(tsNstar,[],1);
	ccc.Neff1d = reshape(tsNeff,[],1);
	ccc.hypothesis1d = reshape(hypothesis,[],1);

	if saveflag == 1
		switch seasonbiasflag
		case 0
			if lower(dataflag(1)) == 'd' & dataflag(end) ~= 'A'
			    save([workdirlocal 'MABclimate2Da_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Da_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			elseif lower(dataflag(1)) == 'b' & dataflag(end) ~= 'A'
				save([workdirlocal 'MABclimate2Db_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Db_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			elseif lower(dataflag(1)) == 'd' & dataflag(end) == 'A'
				save([workdirlocal 'MABclimate2DaSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Da_SODA_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			elseif lower(dataflag(1)) == 'b' & dataflag(end) == 'A'
				save([workdirlocal 'MABclimate2DbSODA_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Db_SODA_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			end %switch
		case 1
			if lower(dataflag(1)) == 'd' & dataflag(end) ~= 'A'
				save([workdirlocal 'MABclimate2Da_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Da_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			elseif lower(dataflag(1)) == 'b' & dataflag(end) ~= 'A'
				save([workdirlocal 'MABclimate2Da_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Db_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			elseif lower(dataflag(1)) == 'd' & dataflag(end) == 'A'
				save([workdirlocal 'MABclimate2DaSODA_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrom','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Da_SODA_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			elseif lower(dataflag(1)) == 'b' & dataflag(end) == 'A'
				save([workdirlocal 'MABclimate2DaSODA_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.mat'],'hydrob','ccc','region','season','yyyy','-v7.3')
				xlswrite([workdirlocal 'MAB_summary_table_2Db_SODA_nobiasDG_r' num2str(searchradius) '_p' num2str(p) '.xlsx'], [ccc.variable1d ccc.region1d ccc.season1d ccc.timeperiod1d ccc.m1d ccc.mlo1d ccc.mhi1d ccc.mloeff1d ccc.mhieff1d ccc.N1d ccc.Nstar1d ccc.Neff1d ccc.hypothesis1d]);
			end %if
		end %switch
	end %if saveflag
end %calcflag
