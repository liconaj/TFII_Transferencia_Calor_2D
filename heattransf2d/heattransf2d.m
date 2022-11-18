classdef heattransf2d
    properties
        TempMesh    % Matriz temperaturas con forma pieza
        k1          % conductividad térmica metal
        k2          % conductividad térmica cámara combustión
        q2          % generacion calor cámara combustión
        h1          % convección térmica ambiente (aire)
        h2          % convección térimca refrigeramiento        
        T1          % temperatura ambiente
        T2          % temperatura refrigeramiento
        SystemMatrix     % Matriz A de AT=b
        SystemVector     % Vector b de AT=b
        SystemTemps      % Vector T de AT=b, resultante de A\b
    end
    properties (GetAccess = private)        
        NodeMesh
        dx
        yp, yn
        xp, xn
    end
    methods
        function obj = heattransf2d(NodeMesh, T, k, h, q)
            obj.NodeMesh = NodeMesh;
            if size(T,2) == 2
                obj.T1 = T(1);
                obj.T2 = T(2);
            else
                obj.T1 = T;
            end
            if size(k,2) == 2
                obj.k1 = k(1);
                obj.k2 = k(2);
            else
                obj.k1 = k;
            end
            if size(h,2) == 2
                obj.h1 = h(1);
                obj.h2 = h(2);
            else
                obj.h1 = h;
            end
            if nargin == 5
                obj.q2 = q;
            end
            obj.dx = obj.NodeMesh.dx;
            obj = solvesystem(obj);
            obj = calcmainnodespos(obj);
        end        
        function showimtemps(obj)
            tempmesh = obj.TempMesh(obj.yp:obj.yn,obj.xp:obj.xn);            
            imagesc(tempmesh)
            colorbar
            daspect([1 1 1])            
            %colormap("jet")
        end
        function obj = solvesystem(obj)
            obj.SystemMatrix = zeros(obj.NodeMesh.numData);
            obj.SystemVector = zeros(obj.NodeMesh.numData, 1);
            ySystem = 0;
            for y = 1:obj.NodeMesh.rows
                for x = 1:obj.NodeMesh.cols
                    nodetype = obj.NodeMesh.Data(y, x);
                    % ignorar nodo si no es metal o cámara
                    % combustion
                    if obj.NodeMesh.DataIndex(y,x) == 0
                        continue
                    end
                    ySystem = ySystem + 1;
                    nodescase = obj.NodeMesh.Data(y-1:y+1,x-1:x+1);
                    [coefMatrix, coefVector] = getcoefs(obj, nodetype, nodescase);
                    xSystem = zeros(1, 5);
                    xSystem(1) = obj.NodeMesh.DataIndex(y, x);
                    xSystem(2) = obj.NodeMesh.DataIndex(y, x-1);
                    xSystem(3) = obj.NodeMesh.DataIndex(y-1, x);
                    xSystem(4) = obj.NodeMesh.DataIndex(y, x+1);
                    xSystem(5) = obj.NodeMesh.DataIndex(y+1, x);
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
        function [Th, consth, nameh] = convecconsts(obj, nametype)
            if contains(nametype, "a")
                Th = obj.T1;
                consth = 2 * obj.h1 * obj.dx / obj.k1;
                nameh = "a";
            elseif contains(nametype, "b") 
                Th = obj.T2;
                consth = 2 * obj.h2 * obj.dx / obj.k1;                
                nameh = "b";
            end
        end
        function [constk, namek] = conducconsts(obj, nametype)
            condm = contains(nametype, "m");
            condg = contains(nametype, "g");            
            if condm && ~condg
                namek = "m";
                constk = 0;
            elseif condg && ~condm
                namek = "g";
                constk = obj.q2 * obj.dx^2 / obj.k2;
            elseif condm && condg
                namek = "g";
                constk = 2 * obj.q2 * obj.dx^2 / obj.k2;
            end
        end
        function [coefMatrix, coefVector] = getcoefs(obj, nodetype, nodescase)
            %{
            b = 2 * obj.q2 * obj.dx / obj.k1;
            currentcase = obj.NodeMesh.Cases(nodetype);
            if contains(currentcase, "a")
                T_convec = obj.T1;
                a = 2 * obj.h1 * obj.dx / obj.k1;                
            elseif contains(currentcase, "b") 
                T_convec = obj.T2;
                a = 2 * obj.h2 * obj.dx / obj.k1;                
            end
            switch currentcase
                case "m"                    
                    coefMatrix = [-4 1 1 1 1];
                    coefVector = 0;                   
                case "g"
                    c = obj.q2 * obj.dx^2 / obj.k2;
                    coefMatrix = [-4 1 1 1 1];
                    coefVector = -c;                    
                case "mmee"
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [0 1 2 1], "e", 1, "side");
                    coefMatrix = [-(4+a) adjcoefs];
                    coefVector = -a * T_convec;
                    %}
                    adjcoefs = shiftcoefs(obj, nodescase, [0 -1/(a+4) -2/(a+4) -1/(a+4)], "e", 1, "side");
                    coefMatrix = [1 adjcoefs];
                    coefVector = a * T_convec / (a + 4);
                case "m3e1"
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [1 1 2 2], "e", 1, "corner");
                    coefMatrix = [-(6+a) adjcoefs];
                    coefVector = -a * T_convec;
                    %}
                    adjcoefs = shiftcoefs(obj, nodescase, [-1/(a+6) -1/(a+6) -2/(a+6) -2/(a+6)], "e", 1, "corner");
                    coefMatrix = [1 adjcoefs];
                    coefVector = a * T_convec / (a + 6);
                case "e3m1"                    
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [1 1 0 0], "m", 1, "corner");
                    coefMatrix = [-(2+a) adjcoefs];
                    coefVector = -a * T_convec;
                    %}
                    adjcoefs = shiftcoefs(obj, nodescase, [-1/(a+2) -1/(a+2) 0 0], "e", 1, "corner");
                    coefMatrix = [1 adjcoefs];
                    coefVector = a * T_convec / (a + 2);
                case "mmgg"
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [2 1 0 1], "m", 1, "side");
                    coefMatrix = [-4 adjcoefs];
                    coefVector = -b;
                    %}
                    adjcoefs = shiftcoefs(obj, nodescase, [-1/2 -1/4 0 -1/4], "m", 1, "side");
                    coefMatrix = [1 adjcoefs];
                    coefVector = b / 4;
                case "m3g1"                    
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [1 1 2 2], "g", 1, "corner");
                    coefMatrix = [-6 adjcoefs];
                    coefVector = -b;
                    %}                    
                    coefMatrix = [1 -1/4 -1/4 -1/4 -1/4];
                    coefVector = b / 4;
                case "mmii" 
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [0 1 2 1], "i", 1, "side");
                    coefMatrix = [-4 adjcoefs];
                    coefVector = 0;
                    %}
                    coefMatrix = [1 -1/4 -2/4 -1/4 0];
                    coefVector = 0;
                case "ggii"
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [0 2 4 2], "i", 1, "side");
                    coefMatrix = [-8 adjcoefs];
                    coefVector = -heatconst;
                    %}
                    coefMatrix = [1 -1/4 -2/4 -1/4 0];
                    coefVector = 0;
                case "i2mg"
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [2 2 0 0], "m", 1, "corner");
                    coefMatrix = [-4 adjcoefs];                    
                    coefVector = -heatconst;
                    %}
                    coefMatrix = [1 -1/4 -2/4 -1/4 0];
                    coefVector = 0;
                case "i2me"
                    %{
                    adjcoefs = shiftcoefs(obj, nodescase, [2 2 0 0], "m", 1, "corner");
                    coefMatrix = [-(1+convecconst) adjcoefs];
                    coefVector = -convecconst * T_convec;
                    %}
                    adjcoefs = shiftcoefs(obj, nodescase, [-1/(a/2+2) -1/(a/2+2) 0 0], "m", 1, "corner");
                    coefMatrix = [1 adjcoefs];
                    coefVector = (a/2 * T_convec)/ (a/2+2);
            end
            %}            
            nametype = obj.NodeMesh.Cases(nodetype);                        
            %disp(nametype)
            if ismember(nametype, ["m", "g"])
                [constk, ~] = conducconsts(obj, nametype);
                adjcoefs = [1 1 1 1];
                coefMatrix = [-4, adjcoefs];
                coefVector = -constk;                
            elseif ismember(nametype,"ggmm")
                [constk, namek] = conducconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [1 1 1 1], namek, "s");                
                coefMatrix = [-4, adjcoefs];
                coefVector = -constk;
            elseif ismember(nametype,"gmmm")                
                adjcoefs = [1 1 1 1];
                coefMatrix = [-4, adjcoefs];
                coefVector = 0;
            elseif ismember(nametype,["aamm", "bbmm"])
                [Th, consth, nameh] = convecconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [0 1 2 1], nameh, "s");
                coefMatrix = [-(consth + 4), adjcoefs];
                coefVector = -consth * Th;
            elseif ismember(nametype,["ammm", "bmmm"])
                [Th, consth, nameh] = convecconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [1 1 2 2], nameh, "c");
                coefMatrix = [-(consth + 6), adjcoefs];
                coefVector = -consth * Th;
            elseif ismember(nametype,["aaam", "bbbm"])
                [~, namek] = conducconsts(obj, nametype);
                [Th, consth, ~] = convecconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [1 1 0 0], namek, "c");
                coefMatrix = [-(consth + 2), adjcoefs];
                coefVector = -consth * Th;
            elseif ismember(nametype,["ggii","iimm"])
                [constk, namek] = conducconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [2 1 0 1], namek, "s");
                coefMatrix = [-4, adjcoefs];                
                coefVector = -constk / 2;
            elseif ismember(nametype,["aiim","biim","aaim","bbim"])
                [~, namek] = conducconsts(obj, nametype);
                [Th, consth, ~] = convecconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [2 2 0 0], namek, "c");
                coefMatrix = [-(consth + 4), adjcoefs];                
                coefVector = -consth * Th;
            elseif ismember(nametype,["iiim", "iiig"])
                [~, namek] = conducconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [1 1 0 0], namek, "c");
                coefMatrix = [-2, adjcoefs];    
                coefVector = 0;                
            elseif ismember(nametype, "giim")
                %[constk, namek] = conducconsts(obj, nametype);
                adjcoefs = shiftcoefs(obj, nodescase, [0 1 2 1], "i", "s");
                coefMatrix = [-4, adjcoefs];
                coefVector = 0;
            end
        end
        function shiftsteps = countshiftsteps(obj, nodescase, nodeseekname, seektype, defnodepos)
            nodeseek = obj.NodeMesh.Types(nodeseekname);
            % adjnodes => [1 2 3 4]
            % "side"  . 2 .   "corner"   1 . 2
            %         1 . 3              . . .
            %         . 4 .              4 . 3
            adjnodes = zeros(1, 4);            
            y = 2; x = 2; % posicion central en matriz currentcase (3x3)
            if seektype == "s"
                adjnodes(1) = nodescase(y, x-1);
                adjnodes(2) = nodescase(y-1, x);
                adjnodes(3) = nodescase(y, x+1);
                adjnodes(4) = nodescase(y+1, x);
            elseif seektype == "c"
                adjnodes(1) = nodescase(y-1, x-1);
                adjnodes(2) = nodescase(y-1, x+1);
                adjnodes(3) = nodescase(y+1, x+1);
                adjnodes(4) = nodescase(y+1, x-1);
            end
            seekpos = find(adjnodes == nodeseek);
            shiftsteps = seekpos - defnodepos;
        end
        function shiftedcoefs = shiftcoefs(obj, nodescase, defcoefs, nodeseekname, seektype, defnodepos)
            if nargin == 5                
                defnodepos = 1;            
            end
            shiftsteps = countshiftsteps(obj, nodescase, nodeseekname, seektype, defnodepos);
            shiftedcoefs = circshift(defcoefs, shiftsteps);
        end       
        function obj = createtempmesh(obj)
            obj.TempMesh = NaN(size(obj.NodeMesh.Data));
            for y = 1:obj.NodeMesh.rows
                for x = 1:obj.NodeMesh.cols
                    systempos = obj.NodeMesh.DataIndex(y,x);
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
    end
end