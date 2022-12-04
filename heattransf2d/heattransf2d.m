classdef heattransf2d
    properties
        TempMesh     % Matriz temperaturas con forma pieza
        SystemMatrix % Matriz A de AT=b
        SystemVector % Vector b de AT=b
        SystemTemps  % Vector T de AT=b, resultante de A\b
        NodeParams   % Diccionario de parametros asociados a un nodo
        HeatConvec   % Diccionario de calor en paredes convección
        TmaxConvec   % Diccionario de temperatura máxima en paredes convección
        PerimConvec  % Diccionario de el perimetro de paredes convección
        MeshIndex
        dx
    end
    properties (GetAccess = private)
        NodeMesh
        yp, yn
        xp, xn
        numeqs
        tprop
        L
    end
    methods
        function obj = heattransf2d(NodeMesh)
            arguments
                NodeMesh nodemesh
            end
            obj.NodeMesh = NodeMesh;
            obj.dx = obj.NodeMesh.dx;
            obj.NodeParams = dictionary();
            obj.TmaxConvec = dictionary();
            obj.tprop = 1;
        end
        function showimtemps(obj,Z)
            arguments
                obj
                Z {mustBeNumeric} = 0
            end
            obj = calcmainnodespos(obj);
            tempmesh = obj.TempMesh(obj.yp:obj.yn,obj.xp:obj.xn);
            [m,n] = size(tempmesh);
            xt = (1:n)*obj.dx;
            yt = (1:m)*obj.dx;
            imagesc(xt,yt,obj.tprop * Z + tempmesh)
            %xt = arrayfun(@num2str,get(gca,'xtick')*obj.dx,'un',0);
            %yt = arrayfun(@num2str,get(gca,'ytick')*obj.dx,'un',0);
            %set(gca,'xticklabel',xt,'yticklabel',yt)
            colorbar
            daspect([1 1 1])
        end
        function showcttemps(obj,Z,levels,filled)
            arguments
                obj
                Z {mustBeNumeric} = 0
                levels = 10
                filled = true
            end
            obj = calcmainnodespos(obj);
            tempmesh = obj.TempMesh(obj.yp:obj.yn,obj.xp:obj.xn);
            tempmesh = tempmesh(end:-1:1,:);
            [m,n] = size(tempmesh);
            xt = (1:n)*obj.dx;
            yt = (1:m)*obj.dx;
            if filled
                contourf(xt,yt,obj.tprop * Z + tempmesh,levels,'LineColor','flat')
            else
                contour(xt,yt,obj.tprop * Z + tempmesh,levels)
            end
            colorbar
            set(gca,'Color','k')
            daspect([1 1 1])
        end
        function obj = solvesystem(obj)
            if isempty(obj.MeshIndex)
                obj = calcmeshindex(obj);
            end
            obj.SystemMatrix = zeros(obj.numeqs);
            obj.SystemVector = zeros(obj.numeqs, 1);
            ySystem = 0;
            tic
            fprintf("Calculando coeficientes... ")
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
            fprintf("Hecho. [%0.fs]\n",toc)
            tic
            fprintf("Resolviendo sistema....... ")
            obj.SystemTemps = obj.SystemMatrix\obj.SystemVector;
            fprintf("Hecho. [%0.fs]\n",toc)
            obj = createtempmesh(obj);
        end
        function obj = setupnk(obj,id,k,q)
            arguments
                obj heattransf2d
                id (1,1) char
                k {mustBeNumeric}
                q {mustBeNumeric} = 0
            end
            ntype = 'K';
            obj.NodeParams(dec2hex(id)) = struct("type",ntype,"k",k,"q",q);
        end
        function obj = setupnh(obj,id,h,T)
            arguments
                obj heattransf2d
                id (1,1) char
                h {mustBeNumeric}
                T {mustBeNumeric}
            end
            ntype = 'H';
            obj.NodeParams(dec2hex(id)) = struct("type",ntype,"h",h,"T",T);
        end
        function obj = setupnq(obj,id,q,A)
            arguments
                obj heattransf2d
                id (1,1) char
                q {mustBeNumeric}
                A {mustBeNumeric,mustBePositive}
            end
            ntype = 'Q';
            obj.NodeParams(dec2hex(id)) = struct("type",ntype,"q",q,"A",A);
        end
        function obj = setupni(obj,id)
            arguments
                obj heattransf2d
                id (1,1) char
            end
            ntype = 'I';
            obj.NodeParams(dec2hex(id)) = struct("type",ntype);
        end
        function obj = setTprop(obj,id,s,L)
            arguments
                obj
                id (1,1) char
                s struct
                L {mustBeNumeric, mustBePositive} = 1
            end
            id = dec2hex(id);
            obj.L = L;
            mf = s.p * s.C;
            obj = calcheatconvec(obj);
            Q = obj.HeatConvec(id);
            obj.tprop =  Q / (mf * s.cp);
        end
        function Tmaxc = getTmaxc(obj,hid,Z)
            arguments
                obj
                hid (1,1) char
                Z {mustBeNumeric} = 0
            end
            Tmaxc = obj.tprop * Z + obj.TmaxConvec(dec2hex(hid));
        end
        function tprop = getTprop(obj)
            tprop = obj.tprop;
        end
        function Tmax = getTmax(obj,Z)
            arguments
                obj
                Z {mustBeNumeric} = 0
            end
            Tmax = max(obj.tprop * Z + obj.TempMesh, [], "all");
        end
        function Q = getHeatConvec(obj,id,opt)
            arguments
                obj
                id = NaN
                opt = ''
            end
            if isnan(id)
                Q = sum(values(obj.HeatConvec)) * obj.L;
                P = sum(values(obj.PerimConvec));
            else
                id = dec2hex(id);
                Q = obj.HeatConvec(id) * obj.L;
                P = obj.PerimConvec(id);
            end
            if strcmp(opt,'A')
                Q = Q / (P * obj.L);
            elseif strcmp(opt,'L')
                Q = Q / obj.L;
            elseif strcmp(opt,'P')
                Q = Q / P;
            end
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
                case 'K2H0Q2I0' %pared con flujo de calor
                    [constq, nameq] = getconstq(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [0 1 2 1], nameq, "s");
                    coefMatrix = [-4, adjcoefs];
                    coefVector = -constq;
                case 'K3H0Q1I0' %esquina interna con flujo de calor
                    [constq, nameq] = getconstq(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [1 1 2 2], nameq, "c");
                    coefMatrix = [-6, adjcoefs];
                    coefVector = -constq;
                case 'K1H0Q1I2' %esquina semiaislada con flujo de calor
                    [~, namek] = getconstk(obj, nodename);
                    [constq, ~] = getconstq(obj, nodename);
                    adjcoefs = shiftcoefs(obj, nodescase, [2 2 0 0], namek, "c");
                    coefMatrix = [-4, adjcoefs];
                    coefVector = -constq;
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
        function shiftedcoefs = shiftcoefs(obj,nodescase, defcoefs, nodeseekname, seektype, defnodepos)
            if nargin == 5
                defnodepos = 1;
            end
            shiftsteps = obj.countshiftsteps(nodescase, nodeseekname, seektype, defnodepos);
            shiftedcoefs = circshift(defcoefs, shiftsteps);
        end
        function [consth, nameh, Th] = getconsth(obj, nodename)
            h = NaN;
            for n = nodename
                nparams = obj.NodeParams(n);
                if nparams.type == 'H'
                    if isnan(h)
                        h = nparams.h;
                    elseif nparams.h ~= h
                        h = h + nparams.h;
                    end
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
        function [constq, nameq] = getconstq(obj, nodename)
            for n = nodename
                nparams = obj.NodeParams(n);
                if nparams.type == 'Q'
                    q = nparams.q;
                    A = nparams.A;
                    nameq = n;
                elseif nparams.type == 'K'
                    k = nparams.k;
                end
            end
            constq = 2 * (q/A) * obj.dx / k;
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
        function obj = calcheatconvec(obj)
            obj.HeatConvec = dictionary();
            obj.PerimConvec = dictionary();
            for y = 1:obj.NodeMesh.rows
                for x = 1:obj.NodeMesh.cols
                    nodename = obj.NodeMesh.Data{y, x};
                    if obj.MeshIndex(y,x) == 0
                        continue
                    end
                    ncase = getcase(obj, nodename);
                    if contains(ncase, ["H1" "H2" "H3"])
                        [~, nameh, Th] = getconsth(obj, nodename);
                        h = obj.NodeParams(nameh).h;
                        T = obj.TempMesh(y, x);
                        Q = h * obj.dx * (T-Th);
                        if ~isConfigured(obj.HeatConvec) || ~isKey(obj.HeatConvec,nameh)
                            obj.PerimConvec(nameh) = 0;
                            obj.HeatConvec(nameh) = 0;
                            obj.TmaxConvec(nameh) = T;
                        end
                        obj.HeatConvec(nameh) = obj.HeatConvec(nameh) + Q;
                        obj.PerimConvec(nameh) = obj.PerimConvec(nameh) + obj.dx;
                        if T > obj.TmaxConvec(nameh)
                            obj.TmaxConvec(nameh) = T;
                        end
                    end
                end
            end
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
        function [heq, eta, efe] = calcfinheq(t,L,Z,k,h)
            % CALCFINHEQ
            % Parámetros
            %   t = grosor
            %   L = longitud
            %   Z = profundidad o ancho
            %   k = coeficiente de conducción de la pieza
            %   h = coeficiente de convección del entorno
            % Retorno
            %   heq = coeficiente de convección equivalente debi a la aleta
            %   eta = eficiencia de la aleta
            %   efe = efectividad de la aleta
            Ak = t * Z; % area de la base
            P = 2 * (t + Z);
            Ac = L * P; % area de la aleta
            nu = sqrt(h * P / (k * Ak));
            eta = tanh(nu * L) /  (nu * L);
            efe = eta * Ac / Ak;
            heq = (h/2) * (1 + eta * Ac/Ak);
        end
        function Pr = calcPr(v,p,cp,k)
            Pr = v * p * cp / k;
        end
        function h = calchint(k,Dh,Pr,Re,Nulam)
            if Re >= 1e4 % turbulento
                Nu = 0.023 * Re^(4/5) * Pr^(1/3);
            else % laminar
                Nu = Nulam;
            end
            h = Nu * k / Dh;
        end
        function h = calchext(k,L,Pr,ReL)
            if ReL >= 5e5 % turbulento
                Nu = 0.037 * ReL^(4/5) * Pr^(1/3);
            else % laminar
                Nu = 0.664 * ReL^(1/2) * Pr^(1/3);
            end
            h = Nu * k / L;
        end
        function Re = calcRe(u,Dh,v)
            Re = u * Dh / v;
        end
        function ReL = calcReL(u,L,v)
            ReL = u * L / v;
        end
        function Dh = calcDh(a,b)
            % Diámetro hidráulico de un ducto rectangular
            Dh = 2 * a * b / (a + b);
        end
    end
end