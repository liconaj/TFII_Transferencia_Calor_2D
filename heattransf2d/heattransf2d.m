classdef heattransf2d
    properties
        TempMesh     % Matriz temperaturas con forma pieza
        SystemMatrix % Matriz A de AT=b
        SystemVector % Vector b de AT=b
        SystemTemps  % Vector T de AT=b, resultante de A\b
        NodeParams
        MeshIndex
    end
    properties (GetAccess = private)        
        NodeMesh
        dx
        yp, yn
        xp, xn
        numeqs
    end
    methods
        function obj = heattransf2d(NodeMesh)
            arguments
                NodeMesh nodemesh
            end
            obj.NodeMesh = NodeMesh;
            obj.dx = obj.NodeMesh.dx; 
        end
        function obj = setnodeparams(obj, nodeparams)
            arguments
                obj heattransf2d
                nodeparams dictionary
            end
            obj.NodeParams = nodeparams;
            obj = calcmeshindex(obj);
        end
        function showimtemps(obj)
            figure
            obj = calcmainnodespos(obj);
            tempmesh = obj.TempMesh(obj.yp:obj.yn,obj.xp:obj.xn);
            imagesc(tempmesh)
            colorbar
            daspect([1 1 1])            
            colormap("hot")
        end
        function obj = solvesystem(obj)
            obj.SystemMatrix = zeros(obj.numeqs);
            obj.SystemVector = zeros(obj.numeqs, 1);
            ySystem = 0;
            for y = 1:obj.NodeMesh.rows
                for x = 1:obj.NodeMesh.cols
                    nodename = obj.NodeMesh.Data{y, x};
                    if obj.MeshIndex(y,x) == 0
                        continue
                    end
                    ySystem = ySystem + 1;
                    nodescase = obj.NodeMesh.Data(y-1:y+1,x-1:x+1);
                    [coefMatrix, coefVector] = getcoefs(obj, nodename, nodescase);
                    xSystem = zeros(1, 5);
                    xSystem(1) = obj.MeshIndex(y, x);
                    xSystem(2) = obj.MeshIndex(y, x-1);
                    xSystem(3) = obj.MeshIndex(y-1, x);
                    xSystem(4) = obj.MeshIndex(y, x+1);
                    xSystem(5) = obj.MeshIndex(y+1, x);
                    for i = 1:5
                        if xSystem(i) == 0
                            continue
                        end
                        obj.SystemMatrix(ySystem, xSystem(i)) = coefMatrix(i);
                        obj.SystemVector(ySystem) = coefVector;
                    end
                end
            end
            obj.SystemTemps = obj.SystemMatrix\obj.SystemVector;            
            obj = createtempmesh(obj);
        end        
    end
    methods (Access = private)
        function [coefMatrix, coefVector] = getcoefs(obj, nodename, nodescase)
            ncase = getcase(obj,nodename);
            switch ncase
                case 'K1H0Q0I0'  % conducción interna
                    [constk, ~] = getconstk(obj, nodename);
                    adjcoefs = [1 1 1 1];
                    coefMatrix = [-4, adjcoefs];
                    coefVector = -constk;
                case 'K2H2Q0I0'  % pared con convección
                    [consth, nameh, Th] = getconsth(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [0 1 2 1], nameh, "s");
                    coefMatrix = [-(consth + 4), adjcoefs];
                    coefVector = -consth * Th;
                case 'K3H1Q0I0' % esquina interna con convección
                    [consth, nameh, Th] = getconsth(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [1 1 2 2], nameh, "c");
                    coefMatrix = [-(consth + 6), adjcoefs];
                    coefVector = -consth * Th;
                case 'K1H3Q0I0' % esquina externa con convección
                    [~, namek] = getconstk(obj, nodename);
                    [consth, ~, Th] = getconsth(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [1 1 0 0], namek, "c");
                    coefMatrix = [-(consth + 2), adjcoefs];
                    coefVector = -consth * Th;
                case 'K2H0Q0I2' % pared aislada
                    [constk, ~] = getconstk(obj, nodename);
                    namei = getnamei(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [0 1 2 1], namei, "s");
                    coefMatrix = [-4, adjcoefs];                
                    coefVector = -constk / 2;
                case 'K1H0Q0I3' % esquina externa aislada
                    [~, namek] = getconstk(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [1 1 0 0], namek, "c");
                    coefMatrix = [-2, adjcoefs];
                    coefVector = 0;
                case 'K1H1Q0I2' %esquina externa semiaislada
                    [~, namek] = getconstk(obj, nodename);
                    [consth, ~, Th] = getconsth(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [2 2 0 0], namek, "c");
                    coefMatrix = [-(consth + 4), adjcoefs];                
                    coefVector = -consth * Th;
                case 'K1H2Q0I1' %esquina externa semiaislada
                    [~, namek] = getconstk(obj, nodename);
                    [consth, ~, Th] = getconsth(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [2 2 0 0], namek, "c");
                    coefMatrix = [-(consth + 4), adjcoefs];                
                    coefVector = -consth * Th;                
            end
        end
        function ncase = getcase(obj,nodename)
            k = 0; h = 0; q = 0; i = 0;
            for n = nodename
                switch obj.NodeParams(n).type
                    case 'K'
                        k = k + 1;
                    case 'H'
                        h = h + 1;
                    case 'Q'
                        q = q + 1;
                    case 'I'
                        i = i + 1;                    
                end
            end
            if k == 4
                k = 1; % frontera dos nodos k diferentes como uno de los dos
            end
            ncase = sprintf('K%dH%dQ%dI%d',k,h,q,i);
        end
        function ignore = ignorenode(obj,ntype)
            utilcount = 0;
            for n = ntype                
                if ~ismember(n,keys(obj.NodeParams))
                    error("Nodo color %d indefinido",hex2dec(n))
                end
                if obj.NodeParams(n).type == 'K'
                    utilcount = utilcount + 1;
                end          
            end
            ignore = utilcount == 0;            
        end
        function [consth, nameh, Th] = getconsth(obj, nodename)
            for n = nodename
                nparams = obj.NodeParams(n);
                if nparams.type == 'H'
                    h = nparams.h;
                    Th = nparams.T;
                    nameh = n;
                elseif nparams.type == 'K'
                    k = nparams.k;
                end
            end
            consth = 2 * h * obj.dx /  k;
        end
        function [constk, namek] = getconstk(obj, nodename)
            k = 0;
            for n = nodename
                nparams = obj.NodeParams(n);
                if nparams.type == 'K' && nparams.k > k
                    k = nparams.k;
                    q = nparams.q;
                    namek = n;
                end
            end
            constk = q * obj.dx^2 / k;
        end
        function namei = getnamei(obj, nodename)
            for n = nodename
                if obj.NodeParams(n).type == 'I'
                    namei = n;
                    break
                end
            end
        end
        function obj = createtempmesh(obj)
            obj.TempMesh = NaN(size(obj.NodeMesh.Data));
            for y = 1:obj.NodeMesh.rows
                for x = 1:obj.NodeMesh.cols
                    systempos = obj.MeshIndex(y,x);
                    if systempos == 0
                        continue
                    end
                    obj.TempMesh(y, x) = obj.SystemTemps(systempos);
                end
            end
        end        
        function obj = calcmainnodespos(obj)
            obj.yp = 1; obj.yn = obj.NodeMesh.rows;
            obj.xp = 1; obj.xn = obj.NodeMesh.cols;
            while all(isnan(obj.TempMesh(:,obj.xp)))
                obj.xp = obj.xp + 1;
            end
            while all(isnan(obj.TempMesh(:,obj.xn)))
                obj.xn = obj.xn - 1;
            end
            while all(isnan(obj.TempMesh(obj.yp,:)))
                obj.yp = obj.yp + 1;
            end
            while all(isnan(obj.TempMesh(obj.yn,:)))
                obj.yn = obj.yn - 1;
            end
        end
        function obj = calcmeshindex(obj)
            [m, n] = size(obj.NodeMesh.Data);
            obj.MeshIndex = zeros(m, n);
            i = 0;
            for y = 1:m
                for x = 1:n
                    ntype = obj.NodeMesh.Data{y,x};
                    if ~ignorenode(obj,ntype)
                        i = i + 1;
                        obj.MeshIndex(y,x) = i;
                    end
                end
            end
            obj.numeqs = i;
        end
        function shiftedcoefs = shiftcoefs(obj,nodescase, defcoefs, nodeseekname, seektype, defnodepos)
            if nargin == 5                
                defnodepos = 1;            
            end
            shiftsteps = obj.countshiftsteps(nodescase, nodeseekname, seektype, defnodepos);
            shiftedcoefs = circshift(defcoefs, shiftsteps);
        end  
    end    
    methods (Static)
        function shiftsteps = countshiftsteps(nodescase, nodeseek, seektype, defnodepos)
            % adjnodes => [1 2 3 4]
            % "side"  . 2 .   "corner"   1 . 2
            %         1 . 3              . . .
            %         . 4 .              4 . 3
            adjnodes = cell(1, 4);            
            y = 2; x = 2; % posicion central en matriz currentcase (3x3)
            if seektype == "s"
                adjnodes{1} = nodescase{y, x-1};
                adjnodes{2} = nodescase{y-1, x};
                adjnodes{3} = nodescase{y, x+1};
                adjnodes{4} = nodescase{y+1, x};
            elseif seektype == "c"
                adjnodes{1} = nodescase{y-1, x-1};
                adjnodes{2} = nodescase{y-1, x+1};
                adjnodes{3} = nodescase{y+1, x+1};
                adjnodes{4} = nodescase{y+1, x-1};
            end
            seekpos = find(string(adjnodes) == nodeseek);
            shiftsteps = seekpos - defnodepos;
        end        
        function nodesetups = setupnk(nodesetups,id,k,q)
            arguments
                nodesetups dictionary
                id (1,1) char
                k {mustBeNumeric,mustBePositive}
                q {mustBeNumeric,mustBeNonnegative} = 0
            end            
            ntype = 'K';
            nodesetups(dec2hex(id)) = struct("type",ntype,"k",k,"q",q);
        end
        function nodesetups = setupnh(nodesetups,id,h,T)
            arguments
                nodesetups dictionary
                id (1,1) char
                h {mustBeNumeric,mustBePositive}
                T {mustBeNumeric,mustBePositive}
            end
            ntype = 'H';
            nodesetups(dec2hex(id)) = struct("type",ntype,"h",h,"T",T);
        end
        function nodesetups = setupnq(nodesetups,id,q,A)
            arguments
                nodesetups dictionary
                id (1,1) char
                q {mustBeNumeric,mustBePositive}
                A {mustBeNumeric,mustBePositive}
            end
            ntype = 'Q';
            nodesetups(dec2hex(id)) = struct("type",ntype,"q",q,"A",A);
        end
        function nodesetups = setupni(nodesetups,id)
            arguments
                nodesetups dictionary
                id (1,1) char            
            end
            ntype = 'I';
            nodesetups(dec2hex(id)) = struct("type",ntype);
        end
    end
end