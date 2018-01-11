function layers = load_splayer(file)  
% Load a picked layer from spicker. 
%
% Developed by DAVE MORSE, University of Washington.

% append our .splayer suffix if it's not already there
% if length(file) > 8
%   suffix = file(length(file)-8:length(file));
%   if ~strcmp(suffix, '.splayer')
%     file=[file '.splayer'];
%   end
% else
%   file=[file '.splayer'];
% end
    
max_nl = 0;
fp=fopen( file, 'rt');
if fp ~= -1
  num=fscanf(fp, '%i', 1);
  for i=1:num,
    np=fscanf(fp, '%i', 1 );
    picks=fscanf(fp, '%i %i %i', [3, np]);
%     xp = picks(1,:);
%     yp = picks(2,:);
%     zp = picks(3,:);
    nl=fscanf(fp, '%i', 1 );
    if nl > max_nl
      max_nl = nl;
    end
    picks=fscanf(fp, '%f %f', [2, nl]);
%     xl = picks(1,:);
%    layers(i,:) = picks(2,:);

  end;
  frewind(fp)
  layers = zeros(num,max_nl).*nan;
  num=fscanf(fp, '%i', 1);
  for i=1:num,
    np=fscanf(fp, '%i', 1 );
    picks=fscanf(fp, '%i %i %i', [3, np]);
%     xp = picks(1,:);
%     yp = picks(2,:);
%     zp = picks(3,:);
    nl=fscanf(fp, '%i', 1 );
    picks=fscanf(fp, '%f %f', [2, nl]);
%     xl = picks(1,:);
   layers(i,1:length(picks)) = picks(2,:);

  end;

  
  fclose( fp );
end;

