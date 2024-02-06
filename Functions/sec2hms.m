function hms = sec2hms(t)
%SEC2HMS Convert seconds into descriptive string.
%   Converts a timespan t in seconds into a string with hours, minutes and 
%   seconds, however not days or even bigger time units. Zero values are 
%   neglected intelligently. The 's' at the end of hour(s), minute(s) and 
%   second(s) is also included dependent on the actual number of time units.
%% Calculate hours, minutes and seconds
hours = floor(t/3600);
t = t - hours * 3600;
mins = floor(t/60);
t = t - mins * 60;
seconds = t;
%% Hours
if hours > 0
    hstr = sprintf('%d hour', hours);
    if hours > 1
        hstr = [hstr 's'];
    end
else
    hstr = '';
end
%% Minutes
if mins > 0
    if hours > 0
        if seconds == 0
            mstr = ' and ';
        else
            mstr = ', ';
        end            
    else
       mstr = '';
    end
	str = sprintf('%d minute', mins);
	mstr = [mstr str];
	if mins > 1
    	mstr = [mstr 's'];
	end
else
    mstr = '';
end
%% Seconds
if seconds > 0
   if mins > 0 || hours > 0
       sstr = ' and ';
   else
       sstr = '';
   end   
   str = sprintf('%.2g second', seconds);
   sstr = [sstr str];
   if seconds ~= 1.0
        sstr = [sstr 's'];
   end
else
    sstr = '';
end
%% Combine all substrings
hms = [hstr mstr sstr];
end