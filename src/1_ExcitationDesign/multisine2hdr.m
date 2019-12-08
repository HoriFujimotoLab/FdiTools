function [path] = multisine2hdr(ms,fname,format)
%MULTISINE2HDR Export reference vector to c header file for multisine with options.
%   [path] = multisine2hdr(ms,fname,format)
% ms        : generated multisine with options
% fname     : header file name e.g. 'data/ms.h'
% format    : output format - 'd' double / 'f' float /
%                             'm' MyWay PE-Expert3 system (MWPE-C6713A DSP)
% Author    : Thomas Beauduin, University of Tokyo
%             Wataru Ohnishi, University of Tokyo
%%%%%
% NOTE:
% Extended from traject2hdr.m in
% https://github.com/HoriFujimotoLab/AcqTools
% Currently, limited to SISO excitation


if nargin < 3
    format = 'd'; % double
end

[folder, name, ext] = fileparts(fname);
fld_nam = [pwd,filesep,folder];
if exist(fld_nam,'dir') == 0, mkdir(fld_nam); end
fid = fopen(strcat(fld_nam,filesep,name,ext),'w');

signal = squeeze(ms.x);
nrofs = size(signal,1);

fprintf(fid, '// CREST FACTOR = %f\n',peak2rms(signal));
fprintf(fid, '// SIGNAL LENGTH = %d [ms]\n\n',nrofs*1e3/ms.harm.fs);
fprintf(fid, '/* HARMONICS PARAMETERS \n');
fprintf(fid, '** ----------------------- \n');
fprintf(fid, '** fs = %f [Hz]: sampling frequency \n', ms.harm.fs);
fprintf(fid, '** df = %f [Hz]: frequency resolution \n', ms.harm.df);
fprintf(fid, '** fl = %f [Hz]: lowest frequency \n', ms.harm.fl);
fprintf(fid, '** fh = %f [Hz]: highest frequency \n', ms.harm.fh);
if strcmp(ms.options.gtp,'q')
    fprintf(fid, '** fr = %f : frequency log ratio \n', ms.harm.fr);
end
fprintf(fid, '*/ \n\n');

fprintf(fid, '/* DESIGN OPTIONS \n');
fprintf(fid, '** ----------------------- \n');
fprintf(fid, '** itp = %s : init phase type:  s=schroeder/r=random \n', ms.options.itp);
fprintf(fid, '** ctp = %s : compression type: c=comp/n=no_comp \n', ms.options.ctp);
fprintf(fid, '** dtp = %s : signal type:      f=full/ O=odd-odd \n**                             o=odd / O2=special odd-odd\n', ms.options.dtp);
fprintf(fid, '** gtp = %s : grid type: l=linear/q=quasi-logarithmic \n', ms.options.gtp);
fprintf(fid, '*/ \n\n');

Hampl = ms.Hampl;
S = whos('Hampl');
if strcmp(S.class,'frd')
    fprintf(fid, '/* NON-PARAMETRIC WEIGHTING */ \n');
else
    [num,den] = tfdata(minreal(ms.Hampl),'v');
    fprintf(fid, '/* AMPLITUDE SPECTRUM \n');
    fprintf(fid, '** ----------------------- \n');
    fprintf(fid, '** tf([num,den]) \n');
    fprintf(fid, '** num = [ ');
    fprintf(fid, '%e ', num);
    fprintf(fid, ']\n');
    fprintf(fid, '** den = [ ');
    fprintf(fid, '%e ', den);
    fprintf(fid, ']\n');
    fprintf(fid, '*/ \n\n');
end

if strcmp(format, 'd') % double
    array_name = ['refvec_' name];
    nrofs_name = ['NROFS_' name];
    fprintf(fid,'#define %s %d \n', nrofs_name, nrofs);
    fprintf(fid,'double %s [%s] = { \n',array_name,nrofs_name);
    
    j = 0;
    for i=1 : nrofs
        fprintf(fid, '%.16e, ', signal(i));
        j=j+1;
        if j == 10
            fprintf(fid, '\n'); j = 0;
        end
    end
    fprintf(fid,'}; \n');
    fclose(fid);
    
elseif strcmp(format, 'f') % float
    array_name = ['refvec_' name];
    nrofs_name = ['NROFS_' name];
    fprintf(fid,'#define %s %d \n', nrofs_name, nrofs);
    fprintf(fid,'float %s [%s] = { \n',array_name,nrofs_name);
    
    j = 0;
    for i=1 : nrofs
        fprintf(fid, '%.8e, ', signal(i));
        j=j+1;
        if j == 10
            fprintf(fid, '\n'); j = 0;
        end
    end
    fprintf(fid,'}; \n');
    fclose(fid);
    
elseif strcmp(format, 'm') % float for myway
    array_name = ['refvec_' name];
    nrofs_name = ['NROFS_' name];
    fprintf(fid,'#define %s %d \n', nrofs_name, nrofs);
    fprintf(fid,'far float %s [%s] = { \n',array_name,nrofs_name);
    
    j = 0;
    for i=1 : nrofs
        fprintf(fid, '%.8e, ', signal(i));
        j=j+1;
        if j == 10
            fprintf(fid, '\n'); j = 0;
        end
    end
    fprintf(fid,'}; \n');
    fclose(fid);
    
else
    error('format: ''d'',''f'',''m''');
end

[folder, name, ext] = fileparts(strcat(fld_nam,filesep,name,ext));
path = strcat(folder, name, ext);

end
