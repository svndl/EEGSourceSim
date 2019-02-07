classdef Node
    % a class for Node in the BrainNetSim
    
    properties
        Type % can be AR, Input , ... non-linear like normalization node
        Freq % ferquency of internal oscillation
        TS % Time seris of a realization for that node
        r %r level defines the pole amplitude, between 0 and 1, closer to 1 more oscillatory
        Number
        Alpha % internal noise level
    end
    
    methods
        
        function this = Node(N)%Type,Freq,TS,rlevel)    
            % class initialization, can be initialized as a vector
            if ~exist('N','var'), N = 1; end
            this(N,1)= this;
            for n = 1:N
                this(n).Type = [];
                this(n).Freq = [];
                this(n).TS = [];
                this(n).r = [];
                this(n).Number = n;
                this(n).Alpha = 1;
            end
        end
    end
    
end