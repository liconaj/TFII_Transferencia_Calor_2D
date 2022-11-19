classdef nodemesh
    properties (GetAccess = public)
        Data          % Matriz de nodos        
        rows          % Número de filas en matriz de nodos
        cols          % Número de columnas en matriz de nodos        
        dx            % [m] distancia dx en malla     
        Sketch

        CellSize        % (m) Tamaño del lado de una celda en sketch (un pixel en la imagen)
        CellDivisions   % Número de partes en la que se va a dividir la celda para crear los nodos   
    end
    properties (GetAccess = private)        
        SketchCols
        SketchRows
        Palette
    end
    properties (Constant)  
        % Paleta de colores de PICO-8
        palettefile = "colormap.png"
    end
    methods
        function obj = nodemesh(file, cellsize, celldivisions)
            arguments
                file string
                cellsize {mustBeNumeric, mustBePositive}
                celldivisions {mustBeInteger, mustBePositive}
            end
            if nargin == 2
                celldivisions = 1;
            end
            obj.CellSize = cellsize;
            obj.CellDivisions = celldivisions;     
            obj.dx = obj.CellSize / obj.CellDivisions;                        
            obj = createsketch(obj, file, celldivisions);
            obj = createnodemesh(obj);           
        end
    end
    methods (Access = private)
        function obj = loadpalette(obj)
            imdata = imread(obj.palettefile);
            [m, n] = size(imdata,[1 2]);
            obj.Palette = zeros(m*n,3);
            i = 0;
            for y = 1:m
                for x = 1:n
                    i = i + 1;
                    color = reshape(imdata(y,x,:),[],3,1);
                    obj.Palette(i,:) = color;
                end
            end
        end
        function obj = createsketch(obj, file, scale)
            obj = loadpalette(obj);
            imdata = imresize(imread(file),scale,"nearest");             
            [m, n] = size(imdata,1:2);
            obj.Sketch = cell(m, n);
            for y = 1:m
                for x = 1:n                    
                    color = reshape(imdata(y,x,:),[],3,1);
                    icolor = find(ismember(obj.Palette, color, "rows"))-1;
                    nodetype = dec2hex(icolor);
                    obj.Sketch{y,x} = nodetype;
                end
            end
            [obj.SketchRows, obj.SketchCols] = size(obj.Sketch, [1 2]);
        end    
        function obj = createnodemesh(obj)
           obj.rows = obj.SketchRows + 1;
           obj.cols = obj.SketchCols + 1;
           obj.Data = cell(obj.rows, obj.cols);             
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
                   typename = obj.gettype(neighbornodes);
                   obj.Data{y, x} = typename;
               end
           end
        end
    end
    methods (Static)
        function typename = gettype(neighbornodes)
            [m, n] = size(neighbornodes);            
            typename = '';
            for i = 1:m
                for j = 1:n
                    ntype = neighbornodes{i,j};
                    typename(end+1) = ntype; %#ok<*AGROW> 
                end
            end
            typename = sort(typename);
            maintype = unique(typename); 
            if length(maintype) == 1
                typename = maintype;
            end
        end
    end
end
