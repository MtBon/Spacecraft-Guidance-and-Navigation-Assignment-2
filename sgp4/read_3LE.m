function [ satrec, longstr1, longstr2 ] = read_3LE( SAT_ID, filename, whichconst )
% Reads TLE number TLENUM in file (or last TLE if TLEUM is larger than number of tle)

    if ~ischar(SAT_ID)
      SAT_ID = num2str(SAT_ID,'%05d');
    end

    fid = fopen(filename);

    % Read line 0
    longstr0 = fgetl(fid);

    while ischar(longstr0)

        % Read first line
        longstr1 = fgetl(fid);

        % Read second line
        longstr2 = fgetl(fid);
        
        
        if strcmp(SAT_ID, longstr1(3:7))

            % Initialize
            typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
            opsmode    = 'a';  % afspc approach
            satrec = twoline2rv( longstr1, longstr2, typerun,'e', opsmode, whichconst);
                        
            break;
        end

        longstr0 = fgetl(fid);

    end

    fclose(fid);

end