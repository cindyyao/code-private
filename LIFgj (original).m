function V_trace = LIFgj(Gexc,Ginh,Gleak,Ggap,Iadd, samplerate,params) 

% Implements a leaky intgrate and fire model, input Gexc and Ginh should be in rows,
% Each row is a cell
% Gleak and Iadd should be a constant, samplerate should be in Hz

% adapted from LIFmodelG: added constant current injection (Iadd)
% JC 5/16/11

% adapted from LIFmodelGplusI: added Ggap, a network of gj conductances
% xyao 2/26/18

% set intrinsic cell parameters (uncomment to run as script)

E_inh = -80 ;%mV      % reversal potential for inhibitory current
E_exc = 0 ;          % reversal potential for excitatory current
E_leak = -60 ;         % reversal potential of leak conductance (this can make cell spike spontaneously)

C = 0.03 ; %nF                          % Set capacitance of cell ref: Diamond 2016 jns DSGC soma size: Sun2002 
V_rest = E_leak ;                          % initial potential of cell (resting potential mV) 
V_threshold = -50 ;                     % spike threshold Graham2001, rabbit DS
abs_ref = .002*samplerate ; %pnts       % 2 ms absolute refractory period (points)
V_relref_tau = .008 ; %sec              % decay time constant of relative refractory period (these numbers taken from Trong and Rieke, 2008)
V_relref_amp = 4 ; %mV                  % amplitude of relative refractory period threshold change 
% Iadd (pA)

% E_inh = params.Einh ;
% E_exc = params.Eexc ;
% E_leak = params.Eleak ;
% C = params.cap ;
% V_rest = params.Vrest ;
% V_threshold = params.Vthresh ;
% abs_ref = round(params.AbsRef * samplerate) ; 
% V_relref_tau = params.RelRefTau ;
% V_relref_amp = params.RelRefAmp ;

% set time (ms)  
time = 1/samplerate:1/samplerate:length(Gexc)/samplerate ;       % time points assessed
tpoints=length(time) ;  % number of time points

V_trace = repmat(time, size(Gexc, 1), 1); 
V_trace(:, 1) = repmat(V_rest, size(Gexc, 1), 1);
I_syn = repmat(time, size(Gexc, 1), 1);
I_syn(:, 1) = zeros(size(Gexc, 1), 1);
I_gap = repmat(time, size(Gexc, 1), 1);
I_gap(:, 1) = zeros(size(Gexc, 1), 1);
ref = zeros(size(Gexc, 1), 1);
V_thr = ones(size(Gexc))*V_threshold;

% trial = 1:size(Gexc,1)
for t = 2:tpoints %for each time point

    for cell = 1:size(Gexc, 1)

        % current from conductance
        I_syn(cell, t) = Gexc(cell, t)*(V_trace(cell,t-1)-E_exc) + Ginh(cell,t)*(V_trace(cell,t-1)-E_inh) + Gleak*(V_trace(cell,t-1)-E_leak) ;
        I_gap(cell, t) = Ggap(cell, :)*(V_trace(cell, t-1) - V_trace(:, t-1));
        if ref(cell) == 0 % if not within refractory period
            V_trace(cell, t) = V_trace(cell, t-1) + (1/samplerate)*((-(I_syn(cell, t) + Iadd + I_gap(cell, t)))/C) ; %h needs to be in units of seconds from ms

        else
            V_trace(cell, t) = V_rest ;
            ref(cell) = ref(cell) - 1 ;

        end

        if V_trace(cell,t)>V_thr(cell,t)
            V_trace(cell,t) = 0 ;
            V_thr(cell,[t:end]) = V_thr(cell,[t:end]) + exp([0:length(V_thr)-t]./(samplerate*-V_relref_tau))*V_relref_amp ; % add on reletive refractory period
            ref(cell) = abs_ref ;
        end

    end % end cell loop

end % end t loop





