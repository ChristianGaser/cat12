function this = circularGraph(adjacencyMatrix,varargin)

      % Constructor
      p = inputParser;
      
      if nargin == 0
        adjacencyMatrix = rand(20);
        thresh = 0.93;
        adjacencyMatrix(adjacencyMatrix >  thresh) = 1;
        adjacencyMatrix(adjacencyMatrix <= thresh) = 0;
        for i = 1:numel(adjacencyMatrix)
          if adjacencyMatrix(i) > 0
            adjacencyMatrix(i) = rand(1,1);
          end
        end
      end
      
      defaultColorMap = jet(length(adjacencyMatrix));
      defaultLabel = cell(length(adjacencyMatrix),1);
      for i = 1:length(defaultLabel)
        defaultLabel{i} = num2str(i);
      end
      
      addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
      
      this.Label    = defaultLabel;
      this.ColorMap = defaultColorMap;
      if nargin > 2
        if strcmp(lower(varargin{1}),'label')
          this.Label    = varargin{2};
        elseif strcmp(lower(varargin{1}),'colormap')
          this.Colormap    = varargin{2};
      end
          
      fig = gcf;
      set(fig,...
        'UserData',this);
      
      try delete(this.Node); end
      clear node;
            
      % Draw the nodes
      t = linspace(-pi,pi,length(adjacencyMatrix) + 1); % theta for each node

      for i = 1:length(adjacencyMatrix)
        this.Node(i) = node(cos(t(i)),sin(t(i)));
        this.Node(i).Color = this.ColorMap(i,:);
        this.Node(i).Label = this.Label{i};
      end

      % Find non-zero values of s and their indices
      [row,col,v] = find(adjacencyMatrix);
      
      % Calculate line widths based on values of s (stored in v).
      minLineWidth  = 0.5;
      lineWidthCoef = 5;
      lineWidth = v./max(v);
      if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
        lineWidth = repmat(minLineWidth,numel(lineWidth),1);
      else % lines of variable width.
        lineWidth = lineWidthCoef*lineWidth + minLineWidth;
      end
      
      % Draw connections on the Poincare hyperbolic disk.
      %
      % Equation of the circles on the disk:
      % x^2 + y^2 
      % + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x 
      % - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
      % where u and v are points on the boundary.
      %
      % Standard form of equation of a circle
      % (x - x0)^2 + (y - y0)^2 = r^2
      %
      % Therefore we can identify
      % x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
      % y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
      % r^2 = x0^2 + y0^2 - 1
      
      for i = 1:length(v)
        if row(i) ~= col(i)
          if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0 
            % points are diametric, so draw a straight line
            u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
            this.Node(row(i)).Connection(end+1) = line(...
              [u(1);v(1)],...
              [u(2);v(2)],...
              'LineWidth', lineWidth(i),...
              'Color', this.ColorMap(row(i),:));
          else % points are not diametric, so draw an arc
            u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0 
              % ensure the arc is within the unit disk
              theta = [linspace(max(thetaLim),pi,50),...
                       linspace(-pi,min(thetaLim),50)];
            else
              theta = linspace(thetaLim(1),thetaLim(2));
            end
            
            this.Node(row(i)).Connection(end+1) = line(...
              r*cos(theta)+x0,...
              r*sin(theta)+y0,...
              'LineWidth', lineWidth(i),...
              'Color', this.ColorMap(row(i),:));
          end
        end
      end
      
      axis image off;
      ax = gca;


      fig = gcf;

      set(fig,'Color',[1 1 1]);
end