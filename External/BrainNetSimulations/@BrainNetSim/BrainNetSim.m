classdef BrainNetSim
    
    properties
        
        Nodes % a strcut array of nodes (NodeNum x 1), storing nodes information including
              % node.Type, Node.Freq, Node.TS %timeSeries
        Connections % a strut matrix of Connections information (NodeNum x NodeNum)including
              % Edge.FiltType, Edge.FiltOrder, Edge.FiltLF,
              % Edge.FiltHF,Edge.Gain
        SF    % sampling frequency for realization of time series
        
        ARMatrix % autoregressive matrix
    end
    
    properties (Dependent)
       NodeNum
    end
    
    methods
        % dependent variable
        function value = get.NodeNum(obj)
            value = numel(obj.Nodes);
        end 
            
        %% network initialization
        function obj = BrainNetSim(NodeNum,SF,rlevel)
            % INPUT: 
                    % NodeNum: number of the nodes in the network
                    % SF: sampling frequency of the system
                    % rLevel: the amplitude of the node poles, value
                    % between 0 and 1, the closer to 1, the more
                    % oscillatory, default value .95
            
            if ~exist('SF','var') || isempty(SF),SF = 100;end
            if ~exist('rlevel','var') || isempty(rlevel),rlevel = 0.95;end
            obj.Nodes = Node(NodeNum);
            [obj.Nodes(:).r] = deal(rlevel);
            obj.Connections = Connection(NodeNum);
            obj.SF = SF;
        end
        
        %% Add a new node
        function [obj,Node_Number] = AddNode (obj,Type,varargin)
            % add an emoty node to the network
            opt = ParseArgs(varargin,...
                'TS'      ,'',...
                'Freq'        ,[],...
                'r'           ,obj.Nodes(end).r, ...
                'Alpha'       , 1 ...
            );
            
            % update node and connection matrices
            Node_Number = obj.NodeNum+1;
            obj.Nodes(Node_Number) = Node(1);
            obj.Connections(Node_Number,Node_Number) = Connection(1);
            
            obj.Nodes(Node_Number).Type = Type;
            obj.Nodes(Node_Number).Freq = opt.Freq;
            obj.Nodes(Node_Number).TS = opt.TS;
            obj.Nodes(Node_Number).r = opt.r;
            obj.Nodes(Node_Number).Number = Node_Number;
            obj.Nodes(Node_Number).Alpha = opt.Alpha;
        end
        
        %% assign internal frequencies
        function obj = AddNodeFreqs (obj,nodes_num, Freqs)
            % add internal frequency to AR nodes
            % syntax NetStrct.AddNodeFreqs(nodes,Freqs)
                % Freqs: is a cell array of internal frequencies of the
                        % nodes, each of the cell should contain one or two
                        % frequencies
            for i = 1:numel(nodes_num)
                % (1) find the node
                ind = find([obj.Nodes(:).Number] == nodes_num(i));
                if ~isempty(Freqs{i})
                    if strcmp(obj.Nodes((ind)).Type,'AR')|| isempty(obj.Nodes((ind)).Type)
                        obj.Nodes((ind)).Type = 'AR';
                        obj.Nodes((ind)).Freq = Freqs{i};
                    else
                        error(['Unable to assign Freq to node' num2str((ind)) '. Check the type of node and change it to AR if necessary.'])
                    end                    
                end
            end
        end
        
        %% add Connections to network
        function obj = AddConnection(obj, nodes,varargin)
            % array of 2x1, indicating the node connections: ATTENTION
            % order is important, it indicates the direction of connection
            % nodes: [source sink] nodes
            opt = ParseArgs(varargin,...
                'Type'      ,'',...
                'LF'        ,[],...
                'HF'        ,[],...
                'Gain'      ,1, ...
                'Order'     ,12 ...
                );
            Connect1 = Connection();
            % check if the connection exist in the network
            if numel(nodes)~=2 || sum(~ismember(nodes, 1:obj.NodeNum))
                error ('nodes should be a 2x1 vector indicating the node number of the network');
            end
            %
            if ~ismember(opt.Type,{'low','high','bandpass','delay'})
                error ('No such filter is defined');
            else
                Connect1.FiltType = opt.Type;
            end
            Connect1.LF = opt.LF;
            Connect1.HF = opt.HF;
            Connect1.Gain = opt.Gain; 
            Connect1.FiltOrder = opt.Order;
            Connect1.Nodes = nodes;
            
            obj.Connections(nodes(1),nodes(2)) = Connect1;
        end
        
        %% AR calculation for linear models only
        function obj = GenerateARMatrix(obj)
            
            % First find the parameters including r(s) and noise level         
            order = [obj.Connections(:).FiltOrder];
            order = max(max(order)+1,4);

            if ~isempty(order) 
            else
                order = 4;
            end
            A = zeros(obj.NodeNum,obj.NodeNum,order);
            
            for n = 1: obj.NodeNum
                f = obj.Nodes(n).Freq/obj.SF; % first frequency
                if ~isempty(f)
                    if numel(f)>1,
                        P = f(2)/f(1);% to have 1/f properties
                    else
                        P = 1;
                    end
                    [r,f] = parameter_search([2*pi*f],obj.Nodes(n).r,P);
                    [x,y] = pol2cart(f,r);Poles = x+1i*y;
                    [~,a] = zp2tf(0,[Poles conj(Poles)],1);
                    A(n,n,1:numel(a)-1)= -a(2:end);
                end
            end

            % the second part is the (filters for connections)
            for n1 = 1:obj.NodeNum
                for n2 = 1:obj.NodeNum
                    if ~isempty(obj.Connections(n1,n2).FiltType) && (n1~=n2) 
                        if sum(ismember(obj.Connections(n1,n2).FiltType,{'high','low','bandpass','stop'}))
                            lf = obj.Connections(n1,n2).LF; hf = obj.Connections(n1,n2).HF;
                            b1 = fir1(obj.Connections(n1,n2).FiltOrder, [lf hf]./(obj.SF/2),obj.Connections(n1,n2).FiltType);
                        elseif strcmp(obj.Connections(n1,n2).FiltType,'delay')
                            A(n1,n2,obj.Connections(n1,n2).FiltOrder) = 1/obj.Connections(n1,n2).Gain;                           
                        else
                            error(['No Filter as ' obj.Connections(n1,n2).FiltType ' is defined']);
                        end
                        A(n1,n2,1:numel(b1))= b1*obj.Connections(n1,n2).Gain;
                    end
                end
            end
            
            obj.ARMatrix = A;
        end
     
        %% Realization
        function [obj,TS] = Realization(obj,NS)
            alpha = 1./(obj.SF/2);   
            order = size(obj.ARMatrix,3);
            if isempty(obj.ARMatrix)
               obj = GenerateARMatrix(obj);
            end
            TS = zeros(obj.NodeNum,order);
            for t = order+1:10000+NS
                TS_temp = zeros(obj.NodeNum,1);
                NAlpha = [obj.Nodes(:).Alpha]';
                for ord = 1:order
                    al = sum(abs(obj.ARMatrix),3);
                    TS_temp = TS_temp + obj.ARMatrix(:,:,ord)'*TS(:,t-ord)+(randn(size(TS_temp,1),1).* NAlpha*alpha);%./(diag(al)+1);
                end
                TS(:,t) = TS_temp;
            end
            TS = TS(:,10001:end);% remove the first samples
            for n = 1:obj.NodeNum
                obj.Nodes(n).TS = TS(n,:);
            end
        end
        
        %% Stability check

    end
    
end

%% Other functions
function [r,f,Gain,F1,P1] = parameter_search(f,r0,P)
% set the amplitude of the poles in a way that gain(f1)/gain(f2)=P
% This needs to be checked: if the number of frequencies are bigger than 2,
% then an optimization algorithm is required

%initial values
f = sort(f); % sort frequencies from low to high
if ~exist('r0','var') || isempty(r0)
    r0 = .95;% = repmat(.95,size(f));
end
r = repmat(r0,size(f));

if numel(f)>1
    for p = numel(f)-1:-1:1
        %  set parameters for search
        P1= 0;
        e = 0.005;% confidence interval
        dis = 1-r(1);flag=0;
        while ((P1>P+e) || (P1<P-e))
            [XP,YP] = pol2cart([f(p:p+1) -f(p:p+1)],[r(p:p+1) r(p:p+1)]);% poles location
            [XT,YT] = pol2cart(f(p:p+1),ones(size(f(p:p+1)))); % the e^jw location
            Dists = squareform(pdist([XP' YP';XT' YT']));
            L = Dists(end-1:end,1:end-2);
            P1 = prod(L(2,:))./prod(L(1,:));
            % search the parameter space
            if flag == 1
                dis = dis/2;
            end
            if P1>P
                r(p) = (r(p)-dis);
            else
                flag=1;
                r(p) = r(p)+dis;
            end
        end
    end
end
F1 = 0:.02:pi;
for F= 1:numel(F1)
    [XP,YP] = pol2cart([f -f],[r r]);% poles location
    [XT,YT] = pol2cart(F1(F),1); % the e^jw location
    Dists = squareform(pdist([XP' YP';XT' YT']));
    L = Dists(end,1:2*numel(f));
    P1(F) = 1./prod(L(1,:));
end
Gain = max(P1);
end


