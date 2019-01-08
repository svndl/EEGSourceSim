%% The Diffusion Model Analysis Toolbox
% Please use the menu on the left to access the DMAT Manual (for the
% command prompt interface), the DMAT GUI Manual (for the graphical
% interface), or the Introduction to Diffusion models. If this is the first
% time you use the DMA Toolbox, please read this document thoroughly (in
% particular the *End User License Agreement*).
%% Introduction
% DMAT is short for the Diffusion Model Analysis Toolbox. It is a MATLAB™
% toolbox for fitting the Ratcliff Diffusion Model to reaction time and
% accuracy data. For a summary of what this means, we refer the user to
% several papers by Roger Ratcliff, in particular Ratcliff (1978,
% Psychological Review, 85, pp.59-108) and Ratcliff and Tuerlinckx (2002,
% Psychonomic Bulletin & Review, 9, pp.438-481).
%%
% A DMAT Primer is to appear in _Behavior Research Methods_ in 2007:
disp(parsetoline(dmatref(3),74))
%
%% Requirements
% To succesfully run DMAT, you need MATLAB version R2006a or better, with
% the Optimization Toolbox and the Statistics Toolbox installed. The
% Symbolic toolbox is nice to have, but DMAT will run without it. Beside
% that, all you need as much RAM and CPU speed as you can get, and perhaps
% a healthy dose of patience and a sense of humor.
%%
% For MATLAB itself, see system requirements at 
% http://www.mathworks.com/support/sysreq/r2006a/.
%
%% Mailing list
% The DMAT mailing list will be used to send information about updates and
% bug fixes, and optionally to handle 'customer service' (although we do
% not guarantee such service). We hope, however, that the mailing list may
% also grow into a forum where researchers who deal with diffusion modeling
% can get in touch with one another easily. We may at some point decide to
% split this subject away from the list for DMAT users, but for now we use
% only one list.
%%
% To subscribe to the mailing list, simply send an e-mail with only the
% words |"subscribe dmatoolbox"| (in the mail body, without quotation
% marks) to listserv@listserv.cc.kuleuven.be. You will receive a
% confirmation mail shortly afterward.
%
%% License summary
% The Diffusion Model Analysis Toolbox (DMAT) is offered free of charge to
% anyone interested in using diffusion models, provided that you properly
% cite it (see below), and that it is not used for financial profit.
% Requests to use DMAT for any commercial purpose must be directed to the
% authors. DMAT comes without any warranty of any kind.
%%
% You are not allowed to redistribute a downloaded copy of DMAT to others.
% If you want others to use DMAT, refer them to
% http://ppw.kuleuven.be/okp/dmatoolbox/.
%% End user license agreement
type dmateula2.txt
%% Citing DMAT
% If you use DMAT for a project that results in a publication, you need to
% properly cite both the software and the accompanying (submitted) article.
% The correct references (in APA style) are:
r1 = dmatref(1);
r2 = dmatref(2);
d1=parsetoline(r1);
d2=parsetoline(r2);
%%
disp(d1)
%%
disp(d2)

%% About the authors
% The DMA Toolbox was created by Joachim Vandekerckhove and Francis
% Tuerlinckx of the Research Group for Quantitative Psychology of the
% Katholieke Universiteit Leuven.

%% Author of this file
% Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
% Part of the DMA Toolbox. Please read the End User License Agreement,
% contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
% See also http://ppw.kuleuven.be/okp/dmatoolbox.