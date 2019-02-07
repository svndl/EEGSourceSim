classdef Connection
    
    properties
        FiltType % high, low, bandpass,stop, delay
        FiltOrder % Order of the filter
        LF % low cut off for lpf and bpf
        HF % high cut off for hpf and bpf
        Gain % Gain of the filter
        Nodes % array of 2x1 indicating input and output Node numbers
    end  
    
    methods
        function this = Connection(N)
            %Class initialization (Empty)
            if ~exist('N','var'), N = 1; end
            this(N,N)= this;
            for n1 = 1:N
                for n2 = 1:N
                    this(n1,n2).FiltType = [];
                    this(n1,n2).FiltOrder = [];
                    this(n1,n2).LF = [];
                    this(n1,n2).HF = [];
                    this(n1,n2).Gain = [];
                    this(n1,n2).Nodes = [];
                end
            end
        end
    end
    
end