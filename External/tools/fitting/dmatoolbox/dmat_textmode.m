function DMAT = dmat_textmode(DMAT)
%DMAT_TEXTMODE  Temporarily leave graphical interface while in DMATGUI
%   DMAT = DMAT_TEXTMODE(DMAT), where DMAT is the main structure in
%   DMATGUI, yields temporary text control.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

fprintf(1,'\n  DMAT is now operating in Text Mode.\n')
fprintf(1,'  Use MATLAB commands as you normally would.\n')
fprintf(1,'  Enter ''textmode=0;'' to return to the graphical interface.\n')
fprintf(1,'  Only changes in the ''DMAT'' structure will be saved.\n\n')
textmode = 1;
while textmode
    eval(input('  DMAT >> ','s'),'beep,disp(lasterr)');
end