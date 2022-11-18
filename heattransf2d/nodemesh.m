classdef nodemesh
    properties (GetAccess = public)
        Data          % Matriz de nodos        
        rows          % Número de filas en matriz de nodos
        cols          % Número de columnas en matriz de nodos
        DataIndex     % Matriz en donde cada elemento es un número que se refiere al índice de un 
                      % elemento en la misma posición en la matriz de nodos que
        numData       % Número de nodos incógnita        
        dx            % [m] distancia dx en malla
        Cases         % Diccionario donde los numeros son llaves y el nombre es el valor
        Types         % Diccionario donde los nombres son llaves y el numero es el valor
        Sketch

        CellSize        % (m) Tamaño del lado de una celda en sketch (un pixel en la imagen)
        CellDivisions   % Número de partes en la que se va a dividir la celda para crear los nodos        
    end
    properties (GetAccess = private)        
        SketchCols
        SketchRows
        MainCases
        AllIgnoreVals
    end
    properties (Constant)
        % I -> Aislamiento (Insulation)
        % A -> Conveccion 1
        % B -> Conveccion 2
        % M -> Metal
        % G -> Generacion
        % F -> Aleta (Fin)

        % NO ALTERAR: el promedio de valores diferentes también suele ser
        % diferente, esto es ideal para que los valores de las fronteras no
        % se repitan
        MainTypes = dictionary(["i" "a" "b" "m" "g" "f"],...
                               [ 11 17  29  43  67  97]);        
        colorlist =...
                [ 0   0   0  ;  % negro
                  0  255 255 ;  % cyan
                  0   0  255 ;  % azul
                 255 255  0  ;  % amarillo
                 255  0   0  ;  % rojo
                  0  255  0  ]; % verde
        IgnoreTypes = ["i", "a", "b"];
    end
    methods
        function obj = nodemesh(file, cellsize, celldivisions)
            if nargin == 2
                celldivisions = 1;
            end            
            if ceil(celldivisions) ~= celldivisions || celldivisions <= 0
                error("parámetro 3 debe ser un entero positivo mayor a 0")
            end            
            obj.CellSize = cellsize;
            obj.CellDivisions = celldivisions;            
            obj.dx = obj.CellSize / obj.CellDivisions;
            obj.MainCases = dictionary(values(obj.MainTypes)', keys(obj.MainTypes)');            
            %fprintf("Valor dx calculado: %0.3f\n", obj.dx)
            obj = createsketch(obj, file, celldivisions);            
            obj.AllIgnoreVals = [];
            for ntype = obj.IgnoreTypes
                obj.AllIgnoreVals(end+1) = obj.MainTypes(ntype);
            end
            obj = createnodemesh(obj);
            obj = createmeshindex(obj);            
        end
    end
    methods (Access = private)
        function obj = createsketch(obj, file, scale)
            imdata = imresize(imread(file),scale,"nearest"); 
            celltypes = values(obj.MainTypes);            
            [m, n] = size(imdata,1:2);
            obj.Sketch = zeros(m, n);
            for i = 1:m
                for j = 1:n
                    r = imdata(i, j, 1);
                    g = imdata(i, j, 2);
                    b = imdata(i, j, 3);
                    color = [r g b];
                    icolor = ismember(obj.colorlist, color, 'rows');
                    obj.Sketch(i,j) = celltypes(icolor);
                end
            end    
            [obj.SketchRows, obj.SketchCols] = size(obj.Sketch, [1 2]);
        end
        function obj = createmeshindex(obj)
            obj.DataIndex = zeros(obj.rows, obj.cols);
            i = 0;
            for y = 1:obj.rows
                for x = 1:obj.cols
                    typeval = obj.Data(y,x);
                    if ~ismember(typeval, obj.AllIgnoreVals)
                        i = i + 1;
                        obj.DataIndex(y, x) = i;                       
                    end
                end
            end
            obj.numData = i;
        end        
        function obj = createnodemesh(obj)
           obj.rows = obj.SketchRows + 1;
           obj.cols = obj.SketchCols + 1;
           obj.Data = zeros(obj.rows, obj.cols);
           typenames = keys(obj.MainTypes)';
           typevalues = values(obj.MainTypes)';    
           for y = 1:obj.rows
               for x = 1:obj.cols
                   switch x
                       case 1
                           xs = x;
                       case obj.cols
                           xs = x - 1;
                       otherwise
                           xs = [x-1, x];
                   end
                   switch y
                       case 1
                           ys = y;
                       case obj.rows
                           ys = y - 1;
                       otherwise
                           ys = [y-1, y];
                   end
                   neighbornodes = obj.Sketch(ys, xs);
                   typeval = mean(neighbornodes, "all");
                   obj.Data(y, x) = typeval;
                   if ~ismember(typeval,typevalues)
                       newtypename = newtype(obj, neighbornodes);
                       typenames(end+1) = newtypename; %#ok<*AGROW> 
                       typevalues(end+1) = typeval;
                       obj = updateignorevals(obj,newtypename, typeval);
                   end
               end
           end
           obj.Cases = dictionary(typevalues,typenames);
           obj.Types = dictionary(typenames,typevalues);
        end
        function obj = updateignorevals(obj, typename, typeval)
            isutil = false;
            for ntype = keys(obj.MainTypes)'
                if ~ismember(ntype, obj.IgnoreTypes) && contains(typename, ntype)
                    isutil = true;
                    break
                end
            end
            if ~isutil
                obj.AllIgnoreVals(end+1) = typeval;
            end
        end
        function typename = newtype(obj,neighbornodes)
            usedtypes = {}; typecount = [];
            [m, n] = size(neighbornodes);
            for i = 1:m
                for j = 1:n
                    nodeval = neighbornodes(i,j);
                    nodename = obj.MainCases(nodeval);
                    usedtypesarray = [usedtypes{:}];
                    if isempty(usedtypesarray) || ~ismember(nodename,usedtypesarray)
                        usedtypes{end+1} = nodename;
                        typecount(end+1) = 1;
                    else
                        postype = find(usedtypesarray==nodename);
                        c = typecount(postype);
                        typecount(postype) = c + 1;
                    end
                end                                
            end
            usedtypes = [usedtypes{:}];
            sortedtypes = sort(usedtypes);
            typename = "";
            if sum(typecount) < 4
                typecount = 2 * typecount;
            end
            for ntype = sortedtypes
                icount = typecount(usedtypes==ntype);
                iname = string(repmat(char(ntype),1,icount));
                typename = sprintf("%s%s",typename,iname);
            end            
        end
        %{
        function caseval = calctypeval(obj, type1,repeats1, type2, repeats2, type3, repeats3)
            sum1 = obj.MainTypes(type1) * repeats1;
            sum2 = obj.MainTypes(type2) * repeats2;
            if nargin == 5
                sum3 = 0;
            else
                sum3 = obj.MainTypes(type3) * repeats3;
            end            
            caseval = (sum1 + sum2 + sum3) / 4;        
        end
        %}        
    end
end
